import os, math, csv
from ij import IJ, ImagePlus, Prefs
from ij import WindowManager as WM
from ij.io import DirectoryChooser, FileSaver
from ij.plugin import HyperStackConverter, ZProjector, ChannelSplitter, Thresholder, Duplicator, RGBStackConverter
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer
from ij.gui import WaitForUserDialog, PointRoi, OvalRoi, NonBlockingGenericDialog
from ij.measure import Measurements

class CellShapeResults(object):
	"""simple class to contain cell shape analysis results"""
	def __init__(self, 
					cell_area_um2=None,
					cell_perimeter_um=None,
					cell_spikiness_index=None, 
					cell_aspect_ratio=None,
					cell_gfp_I_mean=None,
					cell_gfp_I_sd=None, 
					nuclei_in_cell=1
					):
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
	print("Number of cells identified: {}".format(len(rt_areas)));
	for idx in range(len(rt_areas)):
		roim.select(out_imp, idx);
		IJ.run(out_imp, "Set...", "value=255 slice");
	mx_ind = rt_areas.index(max(rt_areas))
	roim.reset();
	roim.close();
	imp.changes = False;
	imp.close();
	return out_imp;


def generate_stats(seg_binary_imp, intensity_channel_imp, cal):
	"""generate output from segmentation image and paired image from which to take intensities"""
	seg_binary_imp.killRoi();
	rt = ResultsTable();
	mxsz = seg_binary_imp.width * seg_binary_imp.height;
	roim = RoiManager();
	pa_options = ParticleAnalyzer.AREA | ParticleAnalyzer.PERIMETER | ParticleAnalyzer.SHAPE_DESCRIPTORS;
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, pa_options, rt, 1000, mxsz);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(seg_binary_imp);
	rt_as = rt.getColumn(rt.getColumnIndex("Area")).tolist();
	rt_ps = rt.getColumn(rt.getColumnIndex("Perim.")).tolist();
	rt_ars = rt.getColumn(rt.getColumnIndex("AR")).tolist();
	I_means = [];
	I_sds = [];
	rois = roim.getRoisAsArray();
	for roi in rois:
		intensity_channel_imp.setRoi(roi);
		stats = intensity_channel_imp.getStatistics(Measurements.MEAN | Measurements.STD_DEV);
		I_means.append(stats.mean);
		I_sds.append(stats.stdDev);
	roim.reset();
	roim.close();
	return [(a * (cal.pixelWidth**2), p * cal.pixelWidth, ar, m, sd) for a, p, ar, m, sd in zip(rt_as, rt_ps, rt_ars, I_means, I_sds)], rois;
	
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
	roim = RoiManager(False);
	for roi in rois:
		roim.addRoi(roi);
	roim.runCommand("Show All with labels")
	RGBStackConverter.convertToRGB(imp);
	roim.moveRoisToOverlay(imp);
	FileSaver(imp).saveAsTiff(output_path);
	roim.runCommand("Show None");
	roim.close();

def save_output_csv(out_statses, output_folder):
	"""save the output of an analysis run to csv"""
	# TODO: deal with this as dictionary!
	f = open(os.path.join(output_folder, "output.csv"), 'wb');
	try:
		writer = csv.writer(f);
		writer.writerow(["Image", "ROI number", "Cell area", "Cell perimeter", "Cell aspect ratio", "GFP channel mean", "GFP channel SD"]);
		for out_stats in out_statses:
			print("out_stats = {}".format(out_stats));
			for roi_idx, stats in enumerate(out_stats[1]):
				print("roi_idx = {}".format(roi_idx));
				print("stast = {}".format(stats));
				writer.writerow([out_stats[0], roi_idx, stats[0], stats[1], stats[2], stats[3], stats[4]]);
	except IOError as e:
		print("problem saving, {}".format(e));
	finally:
		f.close();
	return;

def save_output_csv(out_stats, image_name, output_folder):
	"""generate an output csv on completion of each image rather than at the end of the run"""
	out_path = os.path.join(output_folder, "output.csv");
	csv_exists = os.path.exists(out_path);
	f_open_mode = 'a' if csv_exists else 'wb';
	f = open(out_path, f_open_mode);
	try:
		writer = csv.writer(f);
		if not csv_exists:
			writer.writerow(["Image", "ROI number", "Cell area", "Cell perimeter", "Cell aspect ratio", "GFP channel mean", "GFP channel SD"])
	
		

# SETUP
Prefs.blackBackground = True;
# select folders
dc = DirectoryChooser("choose root folder containing data for analysis");
input_folder = dc.getDirectory();
#dc = DirectoryChooser("choose location to save output");
#output_folder = dc.getDirectory();
#timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
#output_folder = os.path.join(output_folder, (timestamp + ' output'));
#os.mkdir(output_folder);

# load  image(s):
files_lst = [f for f in os.listdir(input_folder) if os.path.splitext(f)[1]=='.tif'];
out_statses = [];
for f in files_lst:
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
	print(metadata);
	cal = imp.getCalibration();
	cal.setUnit(metadata["y_unit"]);
	cal.pixelWidth = metadata["x_physical_size"];
	cal.pixelHeight = metadata["y_physical_size"];
	imp.setCalibration(cal);
	
	# split channels and use gfp to identify which cells have ML1b expression
	channels = ChannelSplitter.split(imp);
	gfp_imp = channels[0];
	gfp_imp.setTitle("GFP");
	threshold_imp = Duplicator().run(gfp_imp);
	threshold_imp.setTitle("GFP_threshold_imp");
	ecad_imp = channels[1];
	ecad_imp.setTitle("E-cadherin");
	nuc_imp = channels[2];
	IJ.run(threshold_imp, "Make Binary", "method=Otsu background=Dark calculate");
	IJ.run(threshold_imp, "Fill Holes", "");
	erode_count = 2;
	for _ in range(erode_count):
		IJ.run(threshold_imp, "Erode", "");
	threshold_imp = keep_blobs_bigger_than(threshold_imp, min_size_pix=1000);
	threshold_imp = my_kill_borders(threshold_imp);
	out_stats, rois = generate_stats(threshold_imp, gfp_imp, cal);
	threshold_imp.changes = False;
	threshold_imp.close();
	print("out_stats = {}".format(out_stats));
	print("len(out_stats) = {}".format(len(out_stats)));
	print("rois = {}".format(rois));
	print("len(rois) = {}".format(len(rois)));
	# save output image
	output_folder = "C:\\Users\\dougk\\Desktop\\dummy output"
	save_qc_image(imp, rois, "{}_plus_overlay.tiff".format(os.path.join(output_folder, os.path.splitext(f)[0])));
	imp.changes = False;
	imp.close();
	out_statses.append((os.path.splitext(f)[0], out_stats));
	# save/append output data
	# get # nuclei per "cell"
save_output_csv(out_statses, output_folder);
	
	