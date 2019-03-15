import os, math
from ij import IJ, ImagePlus
from ij.io import DirectoryChooser
from ij.plugin import HyperStackConverter, ZProjector, ChannelSplitter, Thresholder
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer, GaussianBlur
from ij.gui import WaitForUserDialog, PointRoi, OvalRoi

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

def merge_incorrect_splits_and_get_centroids(imp, centroid_distance_limit=100, size_limit=100):
	"""if particles are found with centroids closer than centroid_distance_limit and both have size<size_limit, get average centroid"""
	imp.killRoi();
	rt = ResultsTable();
	out_imp = IJ.createImage("Nuclei centroids from {}".format(imp.getTitle()), imp.getWidth(), imp.getHeight(), 1, 8);
	out_imp.show();
	IJ.run(out_imp, "Select All", "");
	IJ.run(out_imp, "Set...", "value=0 slice");
	out_imp.show();
	mxsz = imp.width * imp.height;
	roim = RoiManager(False);
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE | ParticleAnalyzer.CENTROID, rt, 0, size_limit);
	pa.setRoiManager(roim);
	roim.reset();
	rt.reset();
	pa.analyze(imp);
#	rt.show("output");
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
	print("new number of nuclei identified = {}".format(len(centroids)));
#	print("centroids = {}".format(centroids));
	roim.reset();
	roim.close();
	for c in centroids:
#		print("c = ({}, {})".format(c[0], c[1]));
		roi = OvalRoi(c[0], c[1], 1, 1);
		out_imp.setRoi(roi);
#		WaitForUserDialog("roi selected...").show();
		IJ.run(out_imp, "Set...", "value=255 slice");
	imp.changes = False;
	imp.close();
	out_imp.show();
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
	imp = HyperStackConverter.toHyperStack(imp, 3, imp.getNSlices()//3, 1, "Color");
	imp = ZProjector.run(imp,"max");
	imp.setC(2);
	IJ.run(imp, "Red", "");
	imp.setC(1);
	IJ.run(imp, "Green", "");
	imp.show();

	# split channels and use ecad to identify which cells have ML1b expression
	channels = ChannelSplitter.split(imp);
	gfp = channels[0];
	ecad = channels[1];
	nuc = channels[2];
	IJ.run(nuc, "Make Binary", "method=Moments background=Dark calculate");
	nuc.show();
#	thr = Thresholder();
#	thr.setBackground("dark");
#	thr.setMethod("IJ_IsoData");
#	nuc_thresh_ip = thr.createMask(nuc);
#	nuc_thresh_ip = nuc.createThresholdMask();
#	nuc_thresh_imp = ImagePlus("Thresholded nuclei", nuc_thresh_ip);
	IJ.run(nuc, "Invert", "");
	IJ.run(nuc, "Fill Holes", "");
	nuc = keep_blobs_bigger_than(nuc, min_size_pix=100);
	IJ.run(nuc, "Watershed", "");
	nuc = keep_blobs_bigger_than(nuc, min_size_pix=1000);
	
	merge_incorrect_splits_and_get_centroids(nuc, centroid_distance_limit=100, size_limit=10000)
#	IJ.run(nuc, "Erode", "stack");
#	nuc_thresh_imp.show()
	

	
	
	
