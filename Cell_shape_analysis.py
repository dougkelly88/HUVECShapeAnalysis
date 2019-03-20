import csv, json, math, os
from datetime import datetime
from ij import IJ, ImagePlus, Prefs
from ij import WindowManager as WM
from ij.io import DirectoryChooser, FileSaver
from ij.plugin import HyperStackConverter, ZProjector, ChannelSplitter, Thresholder, Duplicator, RGBStackConverter
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer
from ij.gui import WaitForUserDialog, PointRoi, OvalRoi, NonBlockingGenericDialog, GenericDialog
from ij.measure import Measurements

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
	_version_string = "0.0.1";

	def __init__(self, last_input_path=None, last_output_path=None, last_analysis_mode=None):
		self.last_input_path = last_input_path;
		self.last_output_path = last_output_path;
		self.last_analysis_mode = last_analysis_mode if last_analysis_mode is not None else self.list_analysis_modes()[0];
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
		return ["GFP intensity", "E-cadherin watershed", "Manual"];

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
					nuclei_in_cell=1
					):
		self.file_name = file_name;
		self.cell_index = cell_index;
		self.cell_area_um2 = cell_area_um2;
		self.cell_perimeter_um = cell_perimeter_um;
		self.cell_spikiness_index = cell_spikiness_index if cell_spikiness_index is not None else self.calculate_cell_spikiness_index(cell_area_um2, cell_perimeter_um);
		self.cell_aspect_ratio = cell_aspect_ratio;
		self.cell_gfp_I_mean = cell_gfp_I_mean;
		self.cell_gfp_I_sd = cell_gfp_I_sd;
		self.nuclei_in_cell = nuclei_in_cell;

	def calculate_cell_spikiness_index(self, cell_area_um2, cell_perimeter_um):
		"""calculate cell spikiness index, i.e. deviation from circular cell"""
		return (cell_perimeter_um**2 / (4 * math.pi * cell_area_um2));

def choose_analysis_mode(params):
	"""present UI for choosing how cells should be identified"""
	dialog = GenericDialog("Analysis methods");
	dialog.addMessage("Please choose how cell shape anlaysis should proceed:");
	dialog.addChoice("Analysis mode", params.list_analysis_modes(), params.last_analysis_mode);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	return dialog.getNextChoice();	

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
	out_imp.show();
	IJ.run(out_imp, "Select All", "");
	IJ.run(out_imp, "Set...", "value=0 slice");
	mxsz = imp.width * imp.height;
	roim = RoiManager();
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, min_size_pix, mxsz);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
	rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
#	print("Number of cells identified: {}".format(len(rt_areas)));
	for idx in range(len(rt_areas)):
		roim.select(out_imp, idx);
		IJ.run(out_imp, "Set...", "value=255 slice");
	mx_ind = rt_areas.index(max(rt_areas))
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

def generate_cell_shape_results(rois, intensity_channel_imp, cal, file_name):
	"""from list of rois, generate results describing the cell enclosed in each roi"""
	pixel_width = 1.0 if cal is None else cal.pixelWidth;
	cell_shapes = [];
	for idx, roi in enumerate(rois):
		intensity_channel_imp.setRoi(roi);
		stats = roi.getStatistics();
		I_mean = stats.mean;
		I_sd = stats.stdDev;
		area = stats.area * (pixel_width**2);
		perimeter = roi.getLength();
		aspect_ratio = stats.major/stats.minor;
		cell_shapes.append(CellShapeResults(file_name=file_name, 
									  cell_index=idx+1, 
									  cell_area_um2=area,
									  cell_perimeter_um=perimeter, 
									  cell_aspect_ratio=aspect_ratio,
									  cell_spikiness_index=None,
									  cell_gfp_I_mean=I_mean,
									  cell_gfp_I_sd=I_sd, 
									  nuclei_in_cell=1));
	return cell_shapes;

def generate_cell_masks(watershed_seeds_imp, intensity_channel_imp):
	"""perform marker-driven watershed on image in intensity_channel_imp"""
	IJ.run(imp, "Marker-controlled Watershed", "input={} marker=Nuclei mask=None binary calculate use".format(os.path.splitext(intensity_channel_imp.getTitle())[0]));
	ws_title =  "{}-watershed".format(intensity_channel_imp.getTitle());
	watershed_imp = WM.getImage(ws_title);
	IJ.setRawThreshold(watershed_imp, 1, watershed_imp.getProcessor().maxValue(), "Red");	
	binary_cells_imp = ImagePlus("thresholded", watershed_imp.createThresholdMask());
	IJ.run(binary_cells_imp, "Kill Borders", "");
	kb_thresh_title = binary_cells_imp.getTitle();
	binary_cells_imp.close();
	binary_cells_imp = WM.getImage("{}-killBorders".format(kb_thresh_title));
	watershed_imp.close();
	return binary_cells_imp;

def merge_incorrect_splits_and_get_centroids(imp, centroid_distance_limit=100, size_limit=100):
	"""if particles are found with centroids closer than centroid_distance_limit and both have size<size_limit, get average centroid"""
	imp.killRoi();
	rt = ResultsTable();
	out_imp = IJ.createImage("Nuclei centroids from {}".format(imp.getTitle()), imp.getWidth(), imp.getHeight(), 1, 8);
	out_imp.show();
	IJ.run(out_imp, "Select All", "");
	IJ.run(out_imp, "Set...", "value=0 slice");
	out_imp.show();
	cal = imp.getCalibration();
	mxsz = imp.width * cal.pixelWidth * imp.height * cal.pixelHeight;
	print("mxsz = {}".format(mxsz));
	roim = RoiManager();
	imp.show();
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE | ParticleAnalyzer.CENTROID, rt, 0, size_limit);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
	MyWaitForUser("paise", "pause post-merge incorrect splits particel analysis");
	rt_xs = rt.getColumn(rt.getColumnIndex("X")).tolist();
	rt_ys = rt.getColumn(rt.getColumnIndex("Y")).tolist();
	centroids = [(x, y) for x, y in zip(rt_xs, rt_ys)];
	print("centroids = {}".format(centroids))
	centroids_set = set();
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
	MyWaitForUser("paise", "pause post-merge incorrect splits particel analysis 2");
	if rt.columnExists("X"):
		rt_xs = rt.getColumn(rt.getColumnIndex("X")).tolist();
		rt_ys = rt.getColumn(rt.getColumnIndex("Y")).tolist();
	centroids = [(x, y) for x, y in zip(rt_xs, rt_ys)];
	for c in centroids:
		centroids_set.add(c);
	centroids = list(centroids_set);
	cal = imp.getCalibration();
	centroids = [(c[0] / cal.pixelWidth, c[1] / cal.pixelHeight) for c in centroids];
	print("new number of nuclei identified = {}".format(len(centroids)));
	roim.reset();
	roim.close();
	for idx, c in enumerate(centroids):
		roi = OvalRoi(c[0], c[1], 10, 10);
		out_imp.setRoi(roi);
		IJ.run(out_imp, "Set...", "value={} slice".format(idx+1));
	imp.changes = False;
	#imp.close();
	return out_imp;

def my_kill_borders(threshold_imp):
	"""handle the tasks required before and after using MorphoLibJ's killBorders"""
	IJ.run(threshold_imp, "Invert", "");
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
					"# nuclei per cell"]);
		for csr in cell_shape_results:
			writer.writerow([csr.file_name, 
								csr.cell_index, 
								csr.cell_area_um2,
								csr.cell_perimeter_um,
								csr.cell_aspect_ratio, 
								csr.cell_spikiness_index, 
								csr.cell_gfp_I_mean, 
								csr.cell_gfp_I_sd, 
								csr.nuclei_in_cell]);
	except IOError as e:
		print("problem saving, {}".format(e));
	finally:
		f.close();
	return;

def gfp_analysis(imp, file_name, output_folder):
	"""perform analysis based on gfp intensity thresholding"""
	cal = imp.getCalibration();
	channel_imps = ChannelSplitter.split(imp);
	gfp_imp = channel_imps[0];
	gfp_imp.setTitle("GFP");
	threshold_imp = Duplicator().run(gfp_imp);
	threshold_imp.setTitle("GFP_threshold_imp");
	ecad_imp = channel_imps[1];
	ecad_imp.setTitle("E-cadherin");
	nuc_imp = channel_imps[2];
	IJ.run(threshold_imp, "Make Binary", "method=Otsu background=Dark calculate");
	IJ.run(threshold_imp, "Fill Holes", "");
	erode_count = 2;
	for _ in range(erode_count):
		IJ.run(threshold_imp, "Erode", "");
	threshold_imp = keep_blobs_bigger_than(threshold_imp, min_size_pix=1000);
	threshold_imp = my_kill_borders(threshold_imp);
	rois = generate_cell_rois(threshold_imp);
	out_stats = generate_cell_shape_results(rois, gfp_imp, cal, file_name);
	print("Number of cells identified = {}".format(len(out_stats)));
	threshold_imp.changes = False;
	threshold_imp.close();
	# save output
	save_qc_image(imp, rois, "{}_plus_overlay.tiff".format(os.path.join(output_folder, os.path.splitext(file_name)[0])));
	save_cell_rois(rois, output_folder, os.path.splitext(file_name)[0])
	imp.changes = False;
	imp.close();
	save_output_csv(out_stats, output_folder);
	return out_stats;

def manual_analysis(imp, file_name, output_folder):
	"""perform analysis based on manually drawn cells"""
	cal = imp.getCalibration();
	channel_imps = ChannelSplitter.split(imp);
	gfp_imp = channel_imps[0];
	IJ.setTool("freehand");
	proceed = False;
	roim = RoiManager();
	roim.runCommand("Show all with labels");
	while not proceed:
		dialog = NonBlockingGenericDialog("Perform manual segmentation");
		dialog.setOKLabel("Proceed to next image...")
		dialog.addMessage("Perform manual segmentation: ");
		dialog.addMessage("Draw around cells and add to the region of interest manager (Ctrl+T)");
		dialog.addMessage("You can see what you've added so far if you check \"show all\" on the ROI manager");
		dialog.addMessage("Then press \"proceed to next image\" when all cells have been added");
		dialog.showDialog();
		if dialog.wasCanceled():
			raise KeyboardInterrupt("Run canceled");
		elif dialog.wasOKed():
			print("rois = {}".format(roim.getCount()));
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
				print("proceeding...")
				rois = roim.getRoisAsArray();
				proceed = True;
	roim.reset();
	roim.close();
	out_stats = generate_cell_shape_results(rois, gfp_imp, cal, file_name);
	print("Number of cells identified = {}".format(len(out_stats)));
	# save output 
	save_qc_image(imp, rois, "{}_plus_overlay.tiff".format(os.path.join(output_folder, os.path.splitext(file_name)[0])));
	save_cell_rois(rois, output_folder, os.path.splitext(file_name)[0])
	imp.changes = False;
	imp.close();
	save_output_csv(out_stats, output_folder);
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
	analysis_mode = choose_analysis_mode(params);
	params.last_analysis_mode = analysis_mode;
	params.persist_parameters();

	# load  image(s):
	files_lst = [f for f in os.listdir(input_folder) if os.path.splitext(f)[1]=='.tif'];
	out_statses = [];
	for f in files_lst:
		print("Working on image {}...".format(os.path.splitext(f)[0]))
		imp = IJ.openImage(os.path.join(input_folder, f));
		metadata = import_iq3_metadata(os.path.join(input_folder, os.path.splitext(f)[0] + '.txt'));
		imp = HyperStackConverter.toHyperStack(imp, 3, imp.getNSlices()//3, 1, "Color");
		imp = ZProjector.run(imp,"max");
		imp.setC(3);
		IJ.run(imp, "Blue", "");
		IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		imp.setC(2);
		IJ.run(imp, "Red", "");
		IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		imp.setC(1);
		IJ.run(imp, "Green", "");
		IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		imp.show();
		imp.setDisplayMode(IJ.COMPOSITE);
		cal = imp.getCalibration();
		cal.setUnit(metadata["y_unit"]);
		cal.pixelWidth = metadata["x_physical_size"];
		cal.pixelHeight = metadata["y_physical_size"];
		imp.setCalibration(cal);
	
		if analysis_mode=="GFP intensity":
			out_stats = gfp_analysis(imp, f, output_folder);
		elif analysis_mode=="Manual":
			out_stats = manual_analysis(imp, f, output_folder);
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

	
	