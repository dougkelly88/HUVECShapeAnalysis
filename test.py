import os, math
from ij import IJ, ImagePlus
from ij import WindowManager as WM
from ij.io import DirectoryChooser
from ij.plugin import HyperStackConverter, ZProjector, ChannelSplitter, Thresholder
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer, GaussianBlur
from ij.gui import WaitForUserDialog, PointRoi, OvalRoi
from ij.measure import Measurements

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
		return meta_dict

def keep_blobs_bigger_than(imp, min_size_pix=100):
	"""remove all blobs other than the largest by area"""
	imp.killRoi();
	rt = ResultsTable();
	if "Size filtered" in imp.getTitle():
		title_addition = "";
	else:
		title_addition = "Size filtered ";
	out_imp = IJ.createImage("{}{}".format(title_addition, imp.getTitle()), imp.getWidth(), imp.getHeight(), 1, 8);
	out_imp.show();
	IJ.run(out_imp, "Select All", "");
	IJ.run(out_imp, "Set...", "value=255 slice");
	mxsz = imp.width * imp.height;
	roim = RoiManager();
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, min_size_pix, mxsz);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
	rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
	print("Number of nuclei identified: {}".format(len(rt_areas)));
	for idx in range(len(rt_areas)):
		roim.select(out_imp, idx);
		IJ.run(out_imp, "Set...", "value=0 slice");
	mx_ind = rt_areas.index(max(rt_areas))
	roim.reset();
	roim.close();
	imp.changes = False;
	imp.close();
	return out_imp;

def generate_stats(seg_binary_imp, intensity_channel_imp, cal):
	"""generate output from segmentation image and paired image from which to take intensities"""
	seg_binary_imp.killRoi();
	IJ.run(seg_binary_imp, "Invert", "");
	rt = ResultsTable();
	mxsz = seg_binary_imp.width * seg_binary_imp.height;
	roim = RoiManager(False);
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
	print(kb_thresh_title);
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
	roim = RoiManager(False);
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE | ParticleAnalyzer.CENTROID, rt, 0, size_limit);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
	rt_xs = rt.getColumn(rt.getColumnIndex("X")).tolist();
	rt_ys = rt.getColumn(rt.getColumnIndex("Y")).tolist();
	centroids = [(x, y) for x, y in zip(rt_xs, rt_ys)];
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
		roi = OvalRoi(c[0], c[1], 1, 1);
		out_imp.setRoi(roi);
		IJ.run(out_imp, "Set...", "value={} slice".format(idx+1));
	imp.changes = False;
	imp.close();
	return out_imp;

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
for f in files_lst:
	imp = IJ.openImage(os.path.join(input_folder, f));
	metadata = import_iq3_metadata(os.path.join(input_folder, os.path.splitext(f)[0] + '.txt'));
	imp = HyperStackConverter.toHyperStack(imp, 3, imp.getNSlices()//3, 1, "Color");
	imp = ZProjector.run(imp,"max");
	imp.setC(2);
	IJ.run(imp, "Red", "");
	imp.setC(1);
	IJ.run(imp, "Green", "");
	imp.show();
	print(metadata);
	cal = imp.getCalibration();
	cal.setUnit(metadata["y_unit"]);
	cal.pixelWidth = metadata["x_physical_size"];
	cal.pixelHeight = metadata["y_physical_size"];
	imp.setCalibration(cal);
	
	# split channels and use ecad to identify which cells have ML1b expression
	channels = ChannelSplitter.split(imp);
	gfp_imp = channels[0];
	ecad_imp = channels[1];
	ecad_imp.setTitle("E-cadherin");
	nuc_imp = channels[2];
	IJ.run(nuc_imp, "Make Binary", "method=Moments background=Dark calculate");
	nuc_imp.show();
	IJ.run(nuc_imp, "Invert", "");
	IJ.run(nuc_imp, "Fill Holes", "");
	nuc_imp = keep_blobs_bigger_than(nuc_imp, min_size_pix=100);
	IJ.run(nuc_imp, "Watershed", "");
	nuc_imp = keep_blobs_bigger_than(nuc_imp, min_size_pix=1000);

	ecad_imp.show();
	ws_seed_imp = merge_incorrect_splits_and_get_centroids(nuc_imp, centroid_distance_limit=100, size_limit=10000);
	binary_cells_imp = generate_cell_masks(ws_seed_imp, ecad_imp) # move stuff below into function...
	ecad_imp.close();
	ws_seed_imp.changes = False;
	ws_seed_imp.close();
	
	# then use particle analyser to generate ROIs for each cell, cycle through ROIs and check whether GFP +ve or not, and generate 
	# shape stats (and save)
	out_stats, rois = generate_stats(binary_cells_imp, gfp_imp, cal);
	print(out_stats);