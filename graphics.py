
def generateVennDiagram(filename = None, script_file = None, plot_file = None):
	"""
		TCGA-2H-A9GF
		area1 = 985,  n=4:, n=5;
		area2 = 618,
		area3 = 1140,
		area4 = 849,
		area5 = 1367,
		n12 = 496,
		n13 = 747,
		n14 = 748,
		n15 = 811,
		n23 = 442,
		n24 = 495,
		n25 = 432,
		n34 = 658,
		n35 = 734,
		n45 = 675,
		n123 = 434,
		n124 = 470,
		n125 = 420,
		n134 = 629,
		n135 = 686,
		n145 = 637,
		n234 = 426,
		n235 = 406,
		n245 = 407,
		n345 = 602,
		n1234 = 420,
		n1235 = 400,
		n1245 = 401,
		n1345 = 580,
		n2345 = 391,
		n12345 = 386,
	"""
	if filename is None:
		filename = "F:\\Combined_Variants\\output\\muse_mutect_somaticsniper_strelka_varscan_Uniquify.vcf"
	if script_file is None:
		script_file = "F:\\Combined_Variants\\output\\rscript.R"
	if plot_file is None:
		plot_file = '"F:/Combined_Variants/output/venn_diagram.tiff"'
	all_caller_names = ['muse', 'mutect', 'somaticsniper', 'strelka', 'varscan']
	combs = list()
	for i in range(1, 6):
		combs += ['_'.join(i) for i in list(itertools.combinations(all_caller_names, i))]
	callers = {k:0 for k in combs}

	with open(filename, 'r') as vcf_file:
		vcf_reader = vcf.Reader(vcf_file)
		for record in vcf_reader:
			category = record.INFO['set']
			category = '_'.join(sorted([i for i in category.split('-') if 'filter' not in i]))
			if category not in callers: callers[category] = 0
			callers[category] += 1

	caller_sets = {'Intersection': callers['Intersection']}
	for key in combs:
		caller_sets[key] = callers['Intersection']
		superset = set(key.split('_'))
		for subset, value in callers.items():
			subset = set(subset.split('_'))
			if superset.issubset(subset):
				caller_sets[key] += value
	r_script = """
		library(VennDiagram)
		venn.plot <- draw.quintuple.venn(
		area1 = {muse},
		area2 = {mutect},
		area3 = {somaticsniper},
		area4 = {strelka},
		area5 = {varscan},
		n12 = {muse_mutect},
		n13 = {muse_somaticsniper},
		n14 = {muse_strelka},
		n15 = {muse_varscan},
		n23 = {mutect_somaticsniper},
		n24 = {mutect_strelka},
		n25 = {mutect_varscan},
		n34 = {somaticsniper_strelka},
		n35 = {somaticsniper_varscan},
		n45 = {strelka_varscan},
		n123 = {muse_mutect_somaticsniper},
		n124 = {muse_mutect_strelka},
		n125 = {muse_mutect_varscan},
		n134 = {muse_somaticsniper_strelka},
		n135 = {muse_somaticsniper_varscan},
		n145 = {muse_strelka_varscan},
		n234 = {mutect_somaticsniper_strelka},
		n235 = {mutect_somaticsniper_varscan},
		n245 = {mutect_strelka_varscan},
		n345 = {somaticsniper_strelka_varscan},
		n1234 = {muse_mutect_somaticsniper_strelka},
		n1235 = {muse_mutect_somaticsniper_varscan},
		n1245 = {muse_mutect_strelka_varscan},
		n1345 = {muse_somaticsniper_strelka_varscan},
		n2345 = {mutect_somaticsniper_strelka_varscan},
		n12345 = {Intersection},
		category = c("MuSE", "MuTect2", "SomaticSniper", "Strelka", "Varscan"),
		fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		cat.cex = 2,
		margin = 0.05,
		cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
		1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
		ind = TRUE
		);

		# Writing to file
		tiff(filename = {filename}, compression = "lzw");
		grid.draw(venn.plot);
		dev.off();
	""".format(**caller_sets, filename = plot_file)
	if script_file:
		with open(script_file, 'w') as outfile:
			outfile.write(r_script)