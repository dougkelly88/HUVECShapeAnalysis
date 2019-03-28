import csv, json, math, os
from datetime import datetime
from ij import IJ, ImagePlus, Prefs
from ij import WindowManager as WM
from ij.io import DirectoryChooser, FileSaver
from ij.plugin import HyperStackConverter, ZProjector, ChannelSplitter, Thresholder, Duplicator, RGBStackConverter
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer
from ij.gui import WaitForUserDialog, PointRoi, OvalRoi, NonBlockingGenericDialog, GenericDialog, PolygonRoi
from ij.measure import Measurements
from ij.process import AutoThresholder

# string definitions
_um = u'\u00B5m';
_degrees = u'\u00B0';
_squared = u'\u00B2';
_totheminusone = u'\u02C9' + u'\u00B9';
_sigma = u'\u03C3'

class Parameters(object):
	"""class to contain parameters for analysis runs"""
	# constants
	_persist_parameters_filename = "IJ_cell_shape_analysis_params.json"
	_persist_parameters_folder = "IJ_cell_shape_analysis";
	_version_string = "0.0.7";

	_analysis_modes = ["GFP intensity", 
						"GFP intensity + manual correction", 
						"E-cadherin watershed", 
						"E-cadherin watershed + manual correction", 
						"Manual"];

	def __init__(self, last_input_path=None, last_output_path=None, last_analysis_mode=None, last_threshold_method='Otsu', last_minimum_cell_area_um2=105):
		self.last_input_path = last_input_path;
		self.last_output_path = last_output_path;
		self.last_analysis_mode = last_analysis_mode if last_analysis_mode is not None else self.list_analysis_modes()[0];
		self.last_threshold_method = last_threshold_method;
		self.last_minimum_cell_area_um2 = last_minimum_cell_area_um2;
		self.__cellshapeparams__ = True;

	def save_parameters_to_json(self, file_path):
		"""save analysis parameters as a json"""
		f = open(file_path, 'w');
		try:
			json.dump(self.__dict__, f, sort_keys=True);
		except Exception as e:
			print(e.message);
			print("error saving {} to {}!".format(self, file_path));
		finally:
			f.close();

	def populate_parameters_from_dict(self, dct):
		"""from dictionary loaded from file, populate parameters"""
		self.last_input_path = dct["last_input_path"];
		self.last_output_path = dct["last_output_path"];
		self.last_analysis_mode = dct["last_analysis_mode"];
		self.last_threshold_method = dct["last_threshold_method"]
		self.last_minimum_cell_area_um2 = dct["last_minimum_cell_area_um2"]

	def load_parameters_from_json(self, file_path):
		"""load parameters from a JSON file"""
		try:
			f = open(file_path, 'r');
			dct = json.loads(f.read());
			if "__cellshapeparams__" in dct:
				self.populate_parameters_from_dict(dct);
			else:
				raise ValueError("JSON file doesn't translate to cell shape analysis parameters")
		except IOError:
			print("IOError reading from JSON file");
			return False;
		except: 
			return False;
		finally:
			f.close();
		return True;

	def load_last_params(self):
		"""load parameters to allow persistence between analysis runs"""
		success = True;
		try:
			temp_path = self.get_peristence_file_location();
			if temp_path:
				temp_params_path = os.path.join(temp_path, Parameters._persist_parameters_filename);
				if os.path.isfile(temp_params_path):
					success = self.load_parameters_from_json(temp_params_path);
				else:
					success = False;
			else:
				success = False;
		except Exception as e:
			print("Warning: Error loading previous settings, reverting to default...");
			raise e;
			return False;
		if not success:
			print("Warning: Error loading previous settings, reverting to default...");
		return success;

	def persist_parameters(self):
		"""save parameters to allow persistence between analysis runs"""
		temp_path = self.get_peristence_file_location();
		if temp_path:
			temp_params_path = os.path.join(temp_path, Parameters._persist_parameters_filename);
			self.save_parameters_to_json(temp_params_path);

	def get_peristence_file_location(self):
		"""get platform-dependent location in which to save persistence data"""
		try:
			if IJ.isWindows():
				# windows
				temp_path = os.path.join(os.getenv('APPDATA'), Parameters._persist_parameters_folder);
			elif IJ.isMacintosh():
				# mac
				temp_path = os.path.join(os.path.expanduser("~"), "Library", Parameters._persist_parameters_folder);
			else:
				print("currently only configured for Mac or Windows - this should be an easy fix for anyone running under Linux...");
				raise NotImplementedError;
			if not os.path.isdir(temp_path):
				os.mkdir(temp_path);
		except Exception as e:
			print("Error: " + e.message);
			return "";
		return temp_path;

	def __str__(self):
		"""return string representation of the Parameters object"""
		return str(self.__dict__);
	
	def list_analysis_modes(self):
		return Parameters._analysis_modes;

class CellShapeResults(object):
	"""simple class to contain cell shape analysis results"""
	def __init__(self, 
					file_name=None,
					cell_index=None,
					cell_area_um2=None,
					cell_perimeter_um=None,
					cell_spikiness_index=None, 
					cell_aspect_ratio=None,
					cell_gfp_I_mean=None,
					cell_gfp_I_sd=None, 
					nuclear_centroids_in_cell=0, 
					nuclei_enclosed_in_cell=0, 
					roi=None):
		self.file_name = file_name;
		self.cell_index = cell_index;
		self.roi = roi;
		self.cell_area_um2 = cell_area_um2;
		self.cell_perimeter_um = cell_perimeter_um;
		self.cell_spikiness_index = cell_spikiness_index;
		self.cell_aspect_ratio = cell_aspect_ratio;
		self.cell_gfp_I_mean = cell_gfp_I_mean;
		self.cell_gfp_I_sd = cell_gfp_I_sd;
		self.nuclear_centroids_in_cell = nuclear_centroids_in_cell;
		self.nuclei_enclosed_in_cell = nuclei_enclosed_in_cell;

	def calculate_cell_spikiness_index(self, cell_area_um2, cell_perimeter_um):
		"""calculate cell spikiness index, i.e. deviation from circular cell"""
		return (cell_perimeter_um**2 / (4 * math.pi * cell_area_um2));

def choose_analysis_mode(params):
	"""present UI for choosing how cells should be identified"""
	dialog = GenericDialog("Analysis methods");
	dialog.addMessage("Please choose how cell shape anlaysis should proceed:");
	dialog.addChoice("Analysis mode: ", params.list_analysis_modes(), params.last_analysis_mode);
	dialog.addChoice("GFP segmentation method: ", AutoThresholder.getMethods(), params.last_threshold_method);
	dialog.addNumericField("Minimum cell area (um" + _squared + "): ", params.last_minimum_cell_area_um2, 0);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	return dialog.getNextChoice(), dialog.getNextChoice(), dialog.getNextNumber();	

def MyWaitForUser(title, message):
	"""non-modal dialog with option to cancel the analysis"""
	dialog = NonBlockingGenericDialog(title);
	dialog.setCancelLabel("Cancel analysis");
	if type(message) is list:
		for line in message:
			dialog.addMessage(line);
	else:
		dialog.addMessage(message);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	return;

def import_iq3_metadata(metadata_path):
	"""import basic image metadata based on the metadata saved by iQ3 software at acquisition time"""
	import re
	x_fmt_str = 'x \: (?P<x_pixels>\d+\.?\d*) \* (?P<x_physical_size>\d+\.?\d*) \: (?P<x_unit>\w+)';
	y_fmt_str = 'y \: (?P<y_pixels>\d+\.?\d*) \* (?P<y_physical_size>\d+\.?\d*) \: (?P<y_unit>\w+)';
	z_fmt_str = '\s*Repeat Z \- (?P<z_extent>[+-]?\d+\.?\d*) (?P<z_unit>\w+) in (?P<z_pixels>\d+\.?\d*) planes \(centre\)';
	t_fmt_str = '\s*Repeat T \- (?P<n_frames>\d+\.?\d*) times \((?P<frame_interval>\d+\.?\d*) (?P<time_unit>(min|sec))\)';
	c_fmt_str = r"\s*Repeat \- Channel \((?P<raw_channels_str>\w.*)\)";
	format_strings = [x_fmt_str, y_fmt_str, z_fmt_str, t_fmt_str, c_fmt_str];
	
	meta_dict = {}
	metadata_file = open(metadata_path, 'r')
	try:
		for line in metadata_file.readlines():
			for fmt_str in format_strings:
				m = re.match(fmt_str, line)
				if (bool(m)):
					meta_dict.update(m.groupdict())
		p_num = re.compile('[+-]?\d+\.?\d*')
		try:
			iteritems = meta_dict.iteritems();
		except:
			iteritems = meta_dict.items();
		for key, value in iteritems:
			if p_num.match(value):
				try:
					meta_dict[key] = float(value)
				except:
					#print("conversion error for key " + key);
					continue;
	finally:
		metadata_file.close();
		if 'raw_channels_str' in meta_dict:
			ch_list = meta_dict['raw_channels_str'].split(",")
			meta_dict['n_channels'] = len(ch_list);
			meta_dict['channel_list'] = ch_list;
		return meta_dict;

def keep_blobs_bigger_than(imp, min_size_pix=100):
	"""remove all blobs other than the largest by area"""
	imp.killRoi();
	rt = ResultsTable();
	if "Size_filtered_" in imp.getTitle():
		title_addition = "";
	else:
		title_addition = "Size_filtered_";
	out_imp = IJ.createImage("{}{}".format(title_addition, imp.getTitle()), imp.getWidth(), imp.getHeight(), 1, 8);
	#out_imp.show();
	IJ.run(out_imp, "Select All", "");
	IJ.run(out_imp, "Set...", "value=0 slice");
	mxsz = imp.width * imp.height;
	roim = RoiManager(False);
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, min_size_pix, mxsz);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
	if roim.getCount()>0:
		rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
	#	print("Number of cells identified: {}".format(len(rt_areas)));
		for idx in range(len(rt_areas)):
			roim.select(out_imp, idx);
			IJ.run(out_imp, "Set...", "value=255 slice");
	roim.reset();
	roim.close();
	imp.changes = False;
	imp.close();
	return out_imp;

def generate_cell_rois(seg_binary_imp):
	"""generate rois from which cell shape information will be gleaned"""
	seg_binary_imp.killRoi();
	mxsz = seg_binary_imp.width * seg_binary_imp.height;
	roim = RoiManager(False);
	pa_options = ParticleAnalyzer.AREA | ParticleAnalyzer.PERIMETER | ParticleAnalyzer.SHAPE_DESCRIPTORS;
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, pa_options, None, 1000, mxsz);
	pa.setRoiManager(roim);
	roim.reset();
	pa.analyze(seg_binary_imp);
	rois = roim.getRoisAsArray();
	roim.reset();
	roim.close();
	return rois;
	
def save_cell_rois(rois, output_folder, filename):
	"""save cell rois to a *.zip file"""
	if len(rois)>0:
		roim = RoiManager(False)
		for roi in rois:
			if roi is not None:
				roim.addRoi(roi);
		roim.runCommand("Save", os.path.join(output_folder, "{} cell rois.zip".format(filename)));
		roim.close();

def generate_cell_shape_results(rois, intensity_channel_imp, cal, file_name, no_nuclei_centroids=None, no_enclosed_nuclei=None):
	"""from list of rois, generate results describing the cell enclosed in each roi"""
	pixel_width = 1.0 if cal is None else cal.pixelWidth;
	if no_nuclei_centroids is None:
		no_nuclei_centroids = [0 for _ in rois];
	if no_enclosed_nuclei is None:
		no_enclosed_nuclei = [0 for _ in rois];
	cell_shapes = [];
	for idx, roi in enumerate(rois):
		intensity_channel_imp.setRoi(roi);
		stats = roi.getStatistics();
		I_mean = stats.mean;
		I_sd = stats.stdDev;
		area = stats.area * (pixel_width**2);
		perimeter = roi.getLength();
		aspect_ratio = stats.major/stats.minor;
		cvh_poly = roi.getConvexHull();
		if cvh_poly is not None:
			convex_hull_roi = PolygonRoi([x for x in cvh_poly.xpoints], [y for y in cvh_poly.ypoints], PolygonRoi.POLYGON);
		else:
			continue;
		print("roi length = {}".format(roi.getLength()));
		print("convex hull roi length = {}".format(convex_hull_roi.getLength()));
		cell_spikiness_index = roi.getLength()/(pixel_width * convex_hull_roi.getLength());
		cell_shapes.append(CellShapeResults(file_name=file_name, 
									  cell_index=idx+1, 
									  cell_area_um2=area,
									  cell_perimeter_um=perimeter, 
									  cell_aspect_ratio=aspect_ratio,
									  cell_spikiness_index=cell_spikiness_index,
									  cell_gfp_I_mean=I_mean,
									  cell_gfp_I_sd=I_sd, 
									  nuclear_centroids_in_cell=no_nuclei_centroids[idx], 
									  nuclei_enclosed_in_cell=no_enclosed_nuclei[idx], 
									  roi=roi));
	return cell_shapes;

def generate_cell_masks(watershed_seeds_imp, intensity_channel_imp, find_edges=False):
	"""perform marker-driven watershed on image in intensity_channel_imp"""
	title = os.path.splitext(intensity_channel_imp.getTitle())[0];
	intensity_channel_imp.show();
	watershed_seeds_imp.show();
	if find_edges:
		IJ.run(intensity_channel_imp, "Find Edges", "");
	IJ.run(intensity_channel_imp, "Marker-controlled Watershed", "input={} marker=Nuclei mask=None binary calculate use".format(title));
	ws_title =  "{}-watershed.tif".format(title);
	watershed_imp = WM.getImage(ws_title);
	IJ.setRawThreshold(watershed_imp, 1, watershed_imp.getProcessor().maxValue(), "Red");	
	binary_cells_imp = ImagePlus("thresholded", watershed_imp.createThresholdMask());
	IJ.run(binary_cells_imp, "Kill Borders", "");
	kb_thresh_title = binary_cells_imp.getTitle();
	binary_cells_imp.close();
	binary_cells_imp = WM.getImage("{}-killBorders".format(kb_thresh_title));
	watershed_imp.close();
	watershed_seeds_imp.changes = False;
	watershed_seeds_imp.close();
	intensity_channel_imp.changes=False;
	intensity_channel_imp.close();
	return binary_cells_imp;

def merge_incorrect_splits_and_get_centroids(imp, centroid_distance_limit=100, size_limit=100):
	"""if particles are found with centroids closer than centroid_distance_limit and both have size<size_limit, get average centroid"""
	imp.killRoi();
	rt = ResultsTable();
	out_imp = IJ.createImage("Nuclei centroids from {}".format(imp.getTitle()), imp.getWidth(), imp.getHeight(), 1, 8);
	IJ.run(out_imp, "Select All", "");
	IJ.run(out_imp, "Set...", "value=0 slice");
	cal = imp.getCalibration();
	mxsz = imp.width * cal.pixelWidth * imp.height * cal.pixelHeight;
	roim = RoiManager(False);
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE | ParticleAnalyzer.CENTROID, rt, 0, size_limit);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
	centroids_set = set();
	if roim.getCount()>0:
		rt_xs = rt.getColumn(rt.getColumnIndex("X")).tolist();
		rt_ys = rt.getColumn(rt.getColumnIndex("Y")).tolist();
		centroids = [(x, y) for x, y in zip(rt_xs, rt_ys)];
		for c in centroids:
			ds = [math.sqrt((c[0] - cx)**2 + (c[1] - cy)**2) for (cx, cy) in centroids];
			close_mask = [d < centroid_distance_limit for d in ds];
			# if no other centroids are within centroid_distance_limit, add this centroid to the output set
			# otherwise, add the average position of this centroid and those within centroid_distance_limit to the output set
			centroids_set.add((sum([msk * b[0] for msk, b in zip(close_mask, centroids)])/sum(close_mask), 
							sum([msk * b[1] for msk, b in zip(close_mask, centroids)])/sum(close_mask)));
	roim.reset();
	rt.reset();
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE | ParticleAnalyzer.CENTROID, rt, size_limit, mxsz);
	pa.setRoiManager(roim);
	pa.analyze(imp);
	if roim.getCount()>0:
		if rt.columnExists("X"):
			rt_xs = rt.getColumn(rt.getColumnIndex("X")).tolist();
			rt_ys = rt.getColumn(rt.getColumnIndex("Y")).tolist();
		centroids = [(x, y) for x, y in zip(rt_xs, rt_ys)];
		for c in centroids:
			centroids_set.add(c);
	centroids = list(centroids_set);
#	cal = imp.getCalibration();
	centroids = [(c[0] / cal.pixelWidth, c[1] / cal.pixelHeight) for c in centroids];
	roim.reset();
	roim.close();
	for idx, c in enumerate(centroids):
		roi = OvalRoi(c[0], c[1], 1, 1);
		out_imp.setRoi(roi);
		IJ.run(out_imp, "Set...", "value={} slice".format(idx+1));
	#imp.changes = False;
	#imp.close();
	return out_imp, centroids;

def my_kill_borders(threshold_imp):
	"""handle the tasks required before and after using MorphoLibJ's killBorders"""
	threshold_imp.killRoi();
#	IJ.run(threshold_imp, "Invert", "");
	IJ.run(threshold_imp, "Kill Borders", "");
	ti_title = threshold_imp.getTitle();
	threshold_imp.changes = False;
	threshold_imp.close();
	threshold_imp = WM.getImage("{}-killBorders".format(ti_title));
	return threshold_imp;

def save_qc_image(imp, rois, output_path):
	"""save rois overlaid on imp to output_path"""
	imp.killRoi();
	if len(rois)>0:
		roim = RoiManager(False);
		for roi in rois:
			roim.addRoi(roi);
		roim.runCommand("Show All with labels")
		RGBStackConverter.convertToRGB(imp);
		roim.moveRoisToOverlay(imp);
		FileSaver(imp).saveAsTiff(output_path);
		roim.runCommand("Show None");
		roim.close();
	else:
		FileSaver(imp).saveAsTiff(output_path);
	return;

def save_output_csv(cell_shape_results, output_folder):
	"""generate an output csv on completion of each image rather than at the end of the run"""
	out_path = os.path.join(output_folder, "output.csv");
	csv_exists = os.path.exists(out_path);
	f_open_mode = 'ab' if csv_exists else 'wb';
	f = open(out_path, f_open_mode);
	try:
		writer = csv.writer(f);
		if not csv_exists:
			writer.writerow(["Image", 
					"ROI number", 
					("Cell area, " + _um + _squared), 
					"Cell perimeter, " + _um, 
					"Cell aspect ratio", 
					"Cell spikiness index", 
					"GFP channel mean", 
					"GFP channel SD", 
					"# nuclear centroids per cell", 
					"# nuclei fully enclosed per cell"]);
		for csr in cell_shape_results:
			writer.writerow([csr.file_name, 
								csr.cell_index, 
								csr.cell_area_um2,
								csr.cell_perimeter_um,
								csr.cell_aspect_ratio, 
								csr.cell_spikiness_index, 
								csr.cell_gfp_I_mean, 
								csr.cell_gfp_I_sd, 
								csr.nuclear_centroids_in_cell, 
								csr.nuclei_enclosed_in_cell]);
	except IOError as e:
		print("problem saving, {}".format(e));
	finally:
		f.close();
	return;

def get_nuclei_locations(nuc_imp, cal, distance_threshold_um=10, size_threshold_um2=50):
	"""get centroids of nuclei. blobs closer than distance_threshold_um and smaller than size_threshold_um2 are merged"""
	size_limit_pix = size_threshold_um2//(cal.pixelWidth**2);
	distance_limit_pix = distance_threshold_um//(cal.pixelWidth);
	IJ.run(nuc_imp, "Make Binary", "method=Moments background=Dark calculate");
	dilate_count = 2;
	for _ in range(dilate_count):
		IJ.run(nuc_imp, "Dilate", "");
	IJ.run(nuc_imp, "Fill Holes", "");
	nuc_imp = keep_blobs_bigger_than(nuc_imp, min_size_pix=math.ceil(float(size_limit_pix)/3));
	nuc_imp.killRoi();
	pre_watershed_nuc_imp = Duplicator().run(nuc_imp);
	IJ.run(nuc_imp, "Watershed", "");
	ws_seed_imp, centroids = merge_incorrect_splits_and_get_centroids(nuc_imp, centroid_distance_limit=distance_limit_pix, size_limit=size_limit_pix);
	full_nuclei_imp = generate_cell_masks(ws_seed_imp, pre_watershed_nuc_imp, find_edges=True);
	return centroids, full_nuclei_imp;

def get_no_nuclei_in_cell(roi, nuclei_centroids):
	"""for a given cell roi and list of nuclei centroids from the image, return how many nuclear centroids lie within the cell"""
	no_nuclei = 0;
	for c in nuclei_centroids:
		if roi.contains(int(c[0]), int(c[1])):
			no_nuclei += 1;
	return no_nuclei;

def get_no_nuclei_fully_enclosed(roi, full_nuclei_imp, overlap_threshold=0.65):
	"""for a given cell roi and ImagePlus with binary nuclei, return how many nuclei lie ENTIRELY within the cell"""
	bbox = roi.getBounds();
	full_nuclei_imp.setRoi(roi);
	cropped_nuc_imp = full_nuclei_imp.crop();
	roi.setLocation(0, 0);
	cropped_nuc_imp.setRoi(roi);
	cropped_nuc_imp.killRoi();
	roim = RoiManager(False);
	mxsz = cropped_nuc_imp.getWidth() * cropped_nuc_imp.getHeight();
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE | ParticleAnalyzer.CENTROID, None, 0, mxsz);
	pa.setRoiManager(roim);
	pa.analyze(cropped_nuc_imp);
	cell_imp = IJ.createImage("Cell binary", 
						cropped_nuc_imp.getWidth(), 
						cropped_nuc_imp.getHeight(), 
						1, 
						8);
	IJ.run(cell_imp, "Select All", "");
	IJ.run(cell_imp, "Set...", "value=0 slice");
	cell_imp.setRoi(roi);
	IJ.run(cell_imp, "Set...", "value=255 slice");
	no_enclosed_nuclei = 0;
	for idx, nuc_roi in enumerate(roim.getRoisAsArray()):
		test_imp = Duplicator().run(cell_imp);
		test_imp.setRoi(nuc_roi);
		IJ.run(test_imp, "Set...", "value=255 slice");
		test_imp.killRoi();
		IJ.run(test_imp, "Create Selection", "");
		IJ.run(test_imp, "Make Inverse", "");
		test_roi = test_imp.getRoi();
		test_roi_stats = test_roi.getStatistics();
		cell_roi_stats = roi.getStatistics();
		nuc_roi_stats = nuc_roi.getStatistics();
		if test_roi_stats.area < (cell_roi_stats.area + (1-overlap_threshold) * nuc_roi_stats.area): # i.e. if more than (100*overlap_threshold)% of nucleus is inside cell...
			no_enclosed_nuclei += 1;
		test_imp.changes = False;
		test_imp.close();
	roi.setLocation(bbox.getX(), bbox.getY());
	cropped_nuc_imp.changes = False;
	cropped_nuc_imp.close();
	cell_imp.changes = False;
	cell_imp.close();
	return no_enclosed_nuclei;

def filter_cells_by_relative_nuclear_area(rois, full_nuclei_imp, relative_nuclear_area_threshold=0.75):
	"""if more than (100*relative_nuclear_area_threshold)% of cell area is made up of nucleus, discard"""
	out_rois = [];
	for roi in rois:
		full_nuclei_imp.setRoi(roi);
		stats = full_nuclei_imp.getStatistics();
		if stats.mean/255 < relative_nuclear_area_threshold:
			out_rois.append(roi);
	return out_rois;


def perform_manual_qc(imp, rois, important_channel=1):
	"""given cell rois generated by automatic methods, allow user to delete/add/redraw as appropriate"""
	for ch in range(imp.getNChannels()):
		imp.setC(ch+1);
		sat_frac = 0.99 if (ch+1)==important_channel else 0.01;
		IJ.run(imp, "Enhance Contrast", "saturated={}".format(sat_frac));

	imp.setC(important_channel);
	IJ.setTool("freehand");
	proceed = False;
	roim = RoiManager();
	roim.runCommand("Show all with labels");
	for roi in rois:
		roim.addRoi(roi);
	while not proceed:
		dialog = NonBlockingGenericDialog("Perform manual segmentation");
		dialog.setOKLabel("Proceed to next image...")
		dialog.addMessage("Perform manual correction of segmentation: ");
		dialog.addMessage("Draw around cells and add to the region of interest manager (Ctrl+T). ");
		dialog.addMessage("Delete and redraw cells as appropriate. ");
		dialog.addMessage("Then press \"proceed to next image\" when all cells have been added. ");
		dialog.showDialog();
		if dialog.wasCanceled():
			raise KeyboardInterrupt("Run canceled");
		elif dialog.wasOKed():
			if roim.getCount()==0:
				rois = [];
				confirm_dialog = GenericDialog("Continue?");
				confirm_dialog.addMessage("No rois selected in this FOV. Are you sure you want to proceed?")
				confirm_dialog.setOKLabel("Yes, proceed");
				confirm_dialog.setCancelLabel("No, not yet");
				confirm_dialog.showDialog();
				if confirm_dialog.wasOKed():
					proceed = True;
			else:
				rois = roim.getRoisAsArray();
				proceed = True;
	roim.reset();
	roim.close();
	for ch in range(imp.getNChannels()):
		imp.setC(ch+1);
		IJ.run(imp, "Enhance Contrast", "saturated={}".format(0.35));
	imp.setC(important_channel);
	return rois;

def ecad_analysis(imp, file_name, output_folder, gfp_channel_number=1, dapi_channel_number=3, red_channel_number=2, do_manual_qc=False):
	"""perform analysis based on marker-driven watershed of junction labelling image"""
	cal = imp.getCalibration();
	channel_imps = ChannelSplitter.split(imp);
	gfp_imp = channel_imps[gfp_channel_number-1];
	gfp_imp.setTitle("GFP");
	ecad_imp = channel_imps[red_channel_number-1];
	ecad_imp.setTitle("E-cadherin");
	threshold_imp = Duplicator().run(ecad_imp);
	threshold_imp.setTitle("Ecad_threshold_imp");
	nuc_imp = channel_imps[dapi_channel_number-1];
	nuclei_locations, full_nuclei_imp, ws_seed_imp = get_nuclei_locations(nuc_imp, cal, distance_threshold_um=10, size_threshold_um2=100);
	full_nuclei_imp.hide();
	binary_cells_imp = generate_cell_masks(ws_seed_imp, ecad_imp, find_edges=False);
	rois = generate_cell_rois(binary_cells_imp);
	if do_manual_qc:
		rois = perform_manual_qc(imp, rois, important_channel=gfp_channel_number);
	no_nuclei_centroids = [get_no_nuclei_in_cell(roi, nuclei_locations) for roi in rois];
	no_enclosed_nuclei = [get_no_nuclei_fully_enclosed(roi, full_nuclei_imp) for roi in rois];
	full_nuclei_imp.changes = False;
	full_nuclei_imp.close();
	out_stats = generate_cell_shape_results(rois, 
										 gfp_imp, 
										 cal, 
										 file_name, 
										 no_nuclei_centroids=no_nuclei_centroids,
										 no_enclosed_nuclei=no_enclosed_nuclei);
	print("Number of cells identified = {}".format(len(out_stats)));
	# save output
	save_qc_image(imp, rois, "{}_plus_overlay.tiff".format(os.path.join(output_folder, os.path.splitext(file_name)[0])));
	save_cell_rois(rois, output_folder, os.path.splitext(file_name)[0])
	imp.changes = False;
	imp.close();
	save_output_csv(out_stats, output_folder);
	return out_stats;

def gfp_analysis(imp, file_name, output_folder, gfp_channel_number=1, dapi_channel_number=3, red_channel_number=2, threshold_method='Otsu', do_manual_qc=False, min_size_pix=1000):
	"""perform analysis based on gfp intensity thresholding"""
	try:
		cal = imp.getCalibration();
		channel_imps = ChannelSplitter.split(imp);
		gfp_imp = channel_imps[gfp_channel_number-1];
		gfp_imp.setTitle("GFP");
		threshold_imp = Duplicator().run(gfp_imp);
		threshold_imp.setTitle("GFP_threshold_imp");
		ecad_imp = channel_imps[red_channel_number-1];
		ecad_imp.setTitle("E-cadherin");
		nuc_imp = channel_imps[dapi_channel_number-1];
		nuclei_locations, full_nuclei_imp = get_nuclei_locations(nuc_imp, cal, distance_threshold_um=10, size_threshold_um2=100);
		full_nuclei_imp.hide();
		IJ.run(threshold_imp, "Make Binary", "method={} background=Dark calculate".format(threshold_method));
		IJ.run(threshold_imp, "Fill Holes", "");
		erode_count = 2;
		for _ in range(erode_count):
			IJ.run(threshold_imp, "Erode", "");
		for _ in range(erode_count):
			IJ.run(threshold_imp, "Dilate", "");
		threshold_imp = keep_blobs_bigger_than(threshold_imp, min_size_pix);
		threshold_imp = my_kill_borders(threshold_imp);
		rois = generate_cell_rois(threshold_imp);
		threshold_imp.changes = False;
		threshold_imp.close();
		rois = filter_cells_by_relative_nuclear_area(rois, full_nuclei_imp, relative_nuclear_area_threshold=0.75);
		if do_manual_qc:
			rois = perform_manual_qc(imp, rois, important_channel=gfp_channel_number);
		no_nuclei_centroids = [get_no_nuclei_in_cell(roi, nuclei_locations) for roi in rois];
		no_enclosed_nuclei = [get_no_nuclei_fully_enclosed(roi, full_nuclei_imp) for roi in rois];
		full_nuclei_imp.changes = False;
		full_nuclei_imp.close();
		out_stats = generate_cell_shape_results(rois, 
											 gfp_imp, 
											 cal, 
											 file_name, 
											 no_nuclei_centroids=no_nuclei_centroids,
											 no_enclosed_nuclei=no_enclosed_nuclei);
		print("Number of cells identified = {}".format(len(out_stats)));
		# save output
		save_qc_image(imp, rois, "{}_plus_overlay.tiff".format(os.path.join(output_folder, os.path.splitext(file_name)[0])));
		save_cell_rois(rois, output_folder, os.path.splitext(file_name)[0])
		imp.changes = False;
		imp.close();
		save_output_csv(out_stats, output_folder);
	except Exception as e:
		print("Ran into a problem analysing {}: {}. Skipping to next cell...".format(file_name, 
																			   e.message));
		out_stats = [];
		pass;
	return out_stats;

def manual_analysis(imp, file_name, output_folder, gfp_channel_number=1, dapi_channel_number=3, red_channel_number=2, important_channel=1):
	"""perform analysis based on manually drawn cells"""
	try:
		cal = imp.getCalibration();
		channel_imps = ChannelSplitter.split(imp);
		gfp_imp = channel_imps[gfp_channel_number-1];
		nuc_imp = channel_imps[dapi_channel_number-1];
		nuclei_locations, full_nuclei_imp = get_nuclei_locations(nuc_imp, cal, distance_threshold_um=10, size_threshold_um2=100);
		full_nuclei_imp.hide();
		rois = perform_manual_qc(imp, [], important_channel=gfp_channel_number);
		no_nuclei_centroids = [get_no_nuclei_in_cell(roi, nuclei_locations) for roi in rois];
		no_enclosed_nuclei = [get_no_nuclei_fully_enclosed(roi, full_nuclei_imp) for roi in rois];
		full_nuclei_imp.changes = False;
		full_nuclei_imp.close();
		out_stats = generate_cell_shape_results(rois, 
											 gfp_imp, 
											 cal, 
											 file_name, 
											 no_nuclei_centroids=no_nuclei_centroids,
											 no_enclosed_nuclei=no_enclosed_nuclei);
		print("Number of cells identified = {}".format(len(out_stats)));
		for ch in range(imp.getNChannels()):
			imp.setC(ch+1);
			IJ.run(imp, "Enhance Contrast", "saturated={}".format(0.35));
		# save output 
		save_qc_image(imp, rois, "{}_plus_overlay.tiff".format(os.path.join(output_folder, os.path.splitext(file_name)[0])));
		save_cell_rois(rois, output_folder, os.path.splitext(file_name)[0])
		imp.changes = False;
		imp.close();
		save_output_csv(out_stats, output_folder);
	except Exception as e:
		print("Ran into a problem analysing {}: {}. Skipping to next cell...".format(file_name, 
																			   e.message));
		out_stats = [];
	return out_stats;

def main():
	# setup
	Prefs.blackBackground = True;
	params = Parameters();
	params.load_last_params();
	# select folders
	if params.last_input_path is not None:
		DirectoryChooser.setDefaultDirectory(params.last_input_path);
	dc = DirectoryChooser("Choose root folder containing data for analysis");
	input_folder = dc.getDirectory();
	params.last_input_path = input_folder;
	if input_folder is None:
		raise KeyboardInterrupt("Run canceled");
	if params.last_output_path is not None:
		DirectoryChooser.setDefaultDirectory(os.path.dirname(params.last_output_path));
	dc = DirectoryChooser("choose location to save output");
	output_folder = dc.getDirectory();
	timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	output_folder = os.path.join(output_folder, (timestamp + ' output'));
	params.last_output_path = output_folder;
	os.mkdir(output_folder);
	analysis_mode, threshold_method, minimum_cell_area_um2 = choose_analysis_mode(params);
	params.last_analysis_mode = analysis_mode;
	do_manual_qc = "+ manual correction" in analysis_mode;
	params.last_threshold_method = threshold_method;
	params.last_minimum_cell_area_um2 = minimum_cell_area_um2;
	params.persist_parameters();

	# load  image(s):
	files_lst = [f for f in os.listdir(input_folder) if os.path.splitext(f)[1]=='.tif'];
	out_statses = [];
	for f in files_lst:
		print("Working on image {}...".format(os.path.splitext(f)[0]))
		imp = IJ.openImage(os.path.join(input_folder, f));
		metadata = import_iq3_metadata(os.path.join(input_folder, os.path.splitext(f)[0] + '.txt'));
		n_channels = int(metadata['n_channels']);
		gfp_channel_number = [("488" in ch) for ch in metadata['channel_list']].index(True) + 1 if any([("488" in ch) for ch in metadata['channel_list']]) else None;
		dapi_channel_number = [("405" in ch) for ch in metadata['channel_list']].index(True) + 1 if any([("405" in ch) for ch in metadata['channel_list']]) else None;
		red_channel_number = [("561" in ch) for ch in metadata['channel_list']].index(True) + 1 if any([("561" in ch) for ch in metadata['channel_list']]) else None;
		imp = HyperStackConverter.toHyperStack(imp, n_channels, imp.getNSlices()//n_channels, 1, "Color");
		imp = ZProjector.run(imp,"max");
		if dapi_channel_number is not None:
			imp.setC(dapi_channel_number);
			IJ.run(imp, "Blue", "");
			IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		else:
			raise NotImplementedError;
		if red_channel_number is not None:
			imp.setC(red_channel_number);
			IJ.run(imp, "Red", "");
			IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		if gfp_channel_number is not None:
			imp.setC(gfp_channel_number);
			IJ.run(imp, "Green", "");
			IJ.run(imp, "Enhance Contrast", "saturated=0.35");
			imp.setC(gfp_channel_number);
		else:
			raise NotImplementedError;
		imp.show();
		imp.setDisplayMode(IJ.COMPOSITE);
		cal = imp.getCalibration();
		cal.setUnit(metadata["y_unit"]);
		cal.pixelWidth = metadata["x_physical_size"];
		cal.pixelHeight = metadata["y_physical_size"];
		imp.setCalibration(cal);

		#threshold_cell_area_um2 = 105; # hardcoded for now, value set for consistency with early runs. TODO: allow user to define during setup
		min_size_pix = minimum_cell_area_um2/(cal.pixelHeight * cal.pixelWidth);
		if "GFP intensity" in analysis_mode:
			out_stats = gfp_analysis(imp, f, output_folder, gfp_channel_number=gfp_channel_number, dapi_channel_number=dapi_channel_number, threshold_method=threshold_method, do_manual_qc=do_manual_qc, min_size_pix=min_size_pix);
		elif analysis_mode=="Manual":
			out_stats = manual_analysis(imp, f, output_folder, gfp_channel_number=gfp_channel_number, dapi_channel_number=dapi_channel_number, important_channel=gfp_channel_number);
		elif "E-cadherin watershed" in analysis_mode:
			#out_stats = ecad_analysis(imp, f, output_folder, gfp_channel_number=gfp_channel_number, dapi_channel_number=dapi_channel_number, do_manual_qc=do_manual_qc);
			imp.close();
			IJ.showMessage("E-cadherin channel thresholding implementation not yet finished!");
			raise NotImplementedError;
		out_statses.extend(out_stats);
		print("Current total number of cells identified: {}".format(len(out_statses)));
		# get # nuclei per "cell"
	params.save_parameters_to_json(os.path.join(output_folder, "parameters used.json"));
	WaitForUserDialog("Done", "Done, having analysed {} cells in total!".format(len(out_statses))).show();
	return;

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main();
