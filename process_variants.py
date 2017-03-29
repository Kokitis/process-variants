import os
import collections
import sys
import configparser
print(sys.version)
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '0_pipeline_files'))
import vcf
import csv
from pprint import pprint

"""
	1. Harmonize VCF fields
	2. Merge caller VCFs
	3. Convert VCFs to MAF
	4. Combine Patient MAFs
"""
if sys.name == 'nt':
	pass
else:
	PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
	
	# SOMATICSEQ_FOLDER/training and SOMATICSEQ/prediction will be created.
	SOMATICSEQ_FOLDER = ""
	OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")

def Terminal(command, show_output = True):
	os.system(command)

def checkdir(path):
	if not os.path.isdir(path):
		os.mkdir(path)
def loadTSV(filename):
	with open(filename, 'r') as file1:
		reader = list(csv.DictReader(file1, delimiter = '\t'))

	return reader

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

class CombineVariants:
	"""Combine separated indel and snv files """
	def __init__(self, sample, options, output_folder, variants = None):
		if variants is None:
			sample_variants = GetVariantList(sample, options)
		else:
			sample_variants = variants
		modified_variants = self._modify_variants(sample_variants, options, output_folder)

		snv_variants, indel_variants = self.splitIndelSnvVariants(modified_variants) #saved in same folder as parent

		output_basename = os.path.join(output_folder, "{0}_vs_{1}.merged".format(sample['NormalID'], sample['SampleID']))
		self.snvs = self.combineVariants(snv_variants, options, output_basename + '.snp.vcf')
		self.indels = self.combineVariants(indel_variants, options, output_basename + '.indel.vcf')

	
	def _modify_merged_vcf(self, filename):
		basename = os.path.splitext(os.path.basename(filename))[0]
		basename = basename + ".modified.vcf"
		output_file = os.path.join(os.path.dirname(filename), basename)

		reader = vcf.Reader(open(filename, 'r'))
		#pprint(reader.infos)
		reader.infos['VAF'] = reader.formats['FREQ']._replace(type='Float')
		#reader.formats['VAF'] = 
		writer = vcf.Writer(open(output_file, 'w'), reader)
	
		for record in reader:
			VAF = self._getVAF(record)
			record.INFO['VAF'] = VAF

			writer.write_record(record)

		return output_file
	
	def _getVAF(self, record):
		""" Use DP4 instead of DP
			Parameters
			----------
				record: 

			Returns
			-------
				result: dict<>
					* 'alleles': The number of alternate alleles.
					* 'reads': The number of reads used to calculate the vaf.
					* 'vaf': The variant allele frequency, in the range 0 - 100.
		"""
		info_keys = record.INFO.keys()
		sample = [i for i in record.samples if i.sample == 'TUMOR'].pop()
		if 'FREQ' in info_keys:
			VAF = float(sample['FREQ'][:-1])
		elif 'DP' in info_keys and 'AD' in info_keys:
			reads = sample['DP']
			alleles = sample['AD'][1]
			VAF = (alleles / reads) * 100
		elif 'AF' in info_keys:
			VAF = sample['AF']
		elif 'DP4' in info_keys:
			reads = sum(sample['DP4'])
			alleles = sum(sample['DP4'][2:])
			VAF = alleles / reads
		else:
			alleles = [i for i in ['A', 'C', 'G', 'T'] if i != record.REF]
			sample_reads = sum([sample[i+'U'][1] for i in (alleles + [record.REF])])
			sample_alleles = sum([sample[i+'U'][1] for i in alleles])
			if sample_reads == 0: sample_vaf = 0
			else:
				sample_vaf = sample_alleles / sample_reads
			VAF = sample_vaf
		return VAF
		"""
		if caller == 'somaticsniper':
			sample_reads = sum(sample['DP4'])
			sample_alleles = sum(sample['DP4'][2:])
			sample_vaf = sample_alleles / sample_reads

		elif caller == 'mutect':
			#print(sample)
			sample_reads = sample['ALT_F1R2'] + sample['ALT_F2R1'] + sample['REF_F1R2'] + sample['REF_F2R1']
			sample_alleles = sample['ALT_F1R2'] + sample['ALT_F2R1']
			sample_vaf = sample['AF']

		elif caller == 'muse':
			sample_reads = sample['DP']
			sample_alleles = sample['AD'][1]
			sample_vaf = sample_alleles / sample_reads

		elif caller == 'strelka-snv':
			alleles = [i for i in ['A', 'C', 'G', 'T'] if i != sample_ref]
			sample_reads = sum([sample[i+'U'][1] for i in (alleles + [sample_ref])])
			sample_alleles = sum([sample[i+'U'][1] for i in alleles])
			if sample_reads == 0: sample_vaf = 0
			else:
				sample_vaf = sample_alleles / sample_reads

		elif 'varscan' in caller:
			sample_reads = sample['DP']
			sample_alleles = sample['AD']
			sample_vaf = float(sample['FREQ'].strip('%')) / 100
		else:
			sample_reads = 0
			sample_alleles = 0
			sample_vaf = 0
			print("WARNING (getSampleVAF): {0} is not a caller!".format(caller))
		"""

		result = {
			'ALLELES': sample_alleles,
			'READS': sample_reads,
			'VAF': float("{0:.5f}".format(sample_vaf)),
		}
		return result

	def _modify_variants(self, vcfs, options, output_folder):
		""" Some caller outputs are inconsistent and need to be modified.
		"""
		#output_folder = "/home/upmc/Documents/Genomic_Analysis/debug_folder/"
		for caller, vcf_filename in vcfs.items():
			basename = os.path.splitext(os.path.basename(vcf_filename))[0] + ".modified.vcf"
			output_file = os.path.join(output_folder, basename)
			with open(vcf_filename, 'r') as vcf_file:
				reader = vcf.Reader(vcf_file)
				if 'varscan' in caller:
					reader = self._modify_varscan_output(reader, output_file)
				vcfs[caller] = self._copy_vcf(reader, output_file)

		return vcfs
	def _copy_vcf(self, reader, output_file):
		with open(output_file, 'w') as file1:
			writer = vcf.Writer(file1, reader)
			for record in reader:
				#print(record.ALT)
				#print(record.ALT[0], type(record.ALT[0]))
				filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
				if filterOut:
					pass
				else:
					writer.write_record(record)
		return output_file
	def _modify_varscan_output(self, reader, filename):
		reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
		reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
		#writer = vcf.Writer(open(filename, 'w'), reader)
		#for record in reader:
		#	writer.write_record(record)
		return reader
	def combineVariants(self, variants, options, output_file):
		""" Uses GATK CombineVariants to merge the calls from each caller into a single file.
			Parameters
			----------
				variants: dict<caller, path>
					A dictionary linkng each caller to its harmonized output.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Returns
			-------
				Output_filename: string
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.merged.vcf
		"""
		gatk = options['Programs']['GATK']
		reference = options['Reference Files']['reference genome']
		#output_folder = '/home/upmc/Documents/Genomic_Analysis/debug_folder/'
		#output_file = "{0}_vs_{1}.merged.vcf".format(sample['NormalID'], sample['SampleID'])
		#output_file = os.path.join(output_folder, output_file)

		order = "mutect,varscan,strelka,muse,somaticsniper" #ordered by VAF confidence
		#order = ",".join([i for i in order.split(',') if i in variants])
		
		#template = "--variant:{0} {1}"
		#templates = ' \\\n'.join([template.format(c, f) for c, f in variants.items()])
		
		commandA = """java -jar "{gatk}" \
			-T CombineVariants \
			-R "{reference}" \
			--variant:muse "{muse}" \
			--variant:mutect "{mutect}" \
			--variant:varscan "{varscan}" \
			--variant:somaticsniper "{ss}" \
			--variant:strelka "{strelka}" \
			-o "{output}" \
			-genotypeMergeOptions PRIORITIZE \
			-priority {rod}"""
		command = commandA.format(
				gatk = 		gatk,
				reference = reference,
				muse = 		variants['muse'],
				mutect = 	variants['mutect'],
				ss = 		variants['somaticsniper'],
				strelka = 	variants['strelka'],
				varscan = 	variants['varscan'],
				rod = 		order,
				output = 	output_file)

		os.system(command)
		return output_file
	def catVariants(self, left, right):
		""" Combines the SNV and Indel files. Assumes both are saved in the same folder. """
		gatk_program = options['Programs']['GATK']
		reference = options['Reference Files']['reference genome']
		l = os.path.splitext(os.path.splitext(left)[0])[0]
		output_file = l + '.cat.vcf'

		command = """java -cp {GATK} org.broadinstitute.gatk.tools.CatVariants \
		    -R {reference}\
		    -V {left} \
		    -V {right} \
		    -out {output}
		""".format(
			GATK = gatk_program,
			reference = reference,
			left = left,
			right = right,
			output = output_file)
		os.system(command)

		return output_file
	def splitIndelSnvVariants(self, variants):
		""" Splits snvs and indels into separate files.
		"""
		snvs = dict()
		indels = dict()
		for caller, filename in variants.items():
			if 'snv' in caller:
				snvs[caller.replace('-snv', '')] = filename
			elif 'indel' in caller:
				indels[caller.replace('-indel', '')] = filename
			else:
				basename = os.path.splitext(filename)[0]
				with open(filename) as vcf_file:
					reader 			= vcf.Reader(vcf_file)
					snv_writer 		= vcf.Writer(open(basename + '.snv.vcf', 'w'), reader)
					indel_writer 	= vcf.Writer(open(basename + '.indel.vcf', 'w'), reader)
					for record in reader:
						if record.is_snp:
							snv_writer.write_record(record)
						elif record.is_indel:
							indel_writer.write_record(record)
				snvs[caller] = basename + '.snv.vcf'
				indels[caller] = basename + '.indel.vcf'
		return snvs, indels


class HarmonizeVCFsOBS:
	def __init__(self, variants, sample, options, **kwargs):
		"""
			Parameters
			----------
				variants: dict<>
					A dictionary mapping a caller to its output file.
				options: dict<>


			Output
			------
				output_folder/harmonized_vcfs/
					Each vcf file will be modified to include any field present in at least one of the raw vcf files.
					These new fields will be empty if no equivilent value is available in the raw vcf file.
				output_folder/merged_vcfs/
					This will contain the output of GATK CombineVariants. Callers are saved in the 
			
		"""
		#training_type = 'RNA-seq'
		self.callers = None #_getMergedRecord has hardcoded values
		self._define_caller_formats()
		
		output_folder = options['Pipeline Options']['processed vcf folder']
		harmonized_folder = os.path.join(output_folder, "harmonized_vcfs")
		merged_folder = os.path.join(output_folder, "merged_vcfs")
		validation_folder = os.path.join(output_folder, "validated_variants")
		
		checkdir(harmonized_folder)
		checkdir(merged_folder)
		checkdir(validation_folder)
		
		print("Output Folder: ", output_folder)
		print("Harmonizing the vcfs...", flush = True)
		variants = self._harmonizeVCFs(sample, variants, output_folder = harmonized_folder)
		
		print("Merging vcfs...", flush = True)
		self.merged_variants = self.merge_vcfs(sample, options, variants, output_folder = merged_folder) #GATK
		
		print("Generating truthset...", flush = True)
		#Truthset is represented as a dict<"chr:{0}|pos:{1} -> validation_status>"
		#truthset = Truthset(sample, options, self.merged_variants, training_type = training_type, **kwargs).truthsetdict
		self.variant_validation_status = collections.defaultdict(int)#truthset#self._get_truthset(truthset)

		print("Validating Variants...", flush = True)
		merged_vcf = self.validate_variants(sample, self.merged_variants, validation_folder)
		
		print("Formatting as a table...")
		merged_table = self.toTable(merged_vcf, validation_folder)

		#_stats_filename = os.path.join(output_folder, "{0}.callset_validation_status.txt".format(training_type))
		#self._count_validated_variants(merged_table, filename = _stats_filename)
	#----------------------------------- Basic ----------------------------------------
	def _define_caller_formats(self):
		muse_format = ["GT", "DP", "AD", "BQ", "SS", 'READS', 'VAF']
		mutect_format = "GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1:READS:VAF".split(':')
		somaticsniper_format = ["GT", "DP", "FDP", "SDP", "SUBDP", "AU", "CU", "GU", "TU","READS", "VAF"]
		strelka_format = "GT:IGT:DP:DP4:BCOUNT:GQ:JGQ:VAQ:BQ:MQ:AMQ:SS:SSC:READS:VAF".split(":")
		varscan_format = ['GT', 'GQ', 'DP', 'RD', 'AD', 'FREQ', 'DP4', 'READS', 'VAF'] 
		all_formats = muse_format + mutect_format + somaticsniper_format + strelka_format + varscan_format
		all_formats = ['GT'] + sorted([ i for i in set(all_formats) if i != 'GT']) #GT has to be first
		all_formats.remove('AF')
		all_formats.remove('FREQ')
		self.muse_tuple = collections.namedtuple('museData', muse_format)
		self.mutect_tuple = collections.namedtuple('mutectData', mutect_format)
		self.strelka_tuple = collections.namedtuple('strelkaData', strelka_format)
		self.ss_tuple = collections.namedtuple('SSData', somaticsniper_format)
		self.varscan_tuple = collections.namedtuple('varscanData', varscan_format)
		
		self.harmonized_tuple = collections.namedtuple('Format', all_formats)
	
	def _get_merged_data(self, fields):
		""" """
		harm_fields = list(self.harmonized_tuple._fields)
		for f in list(fields.keys()):
			if f not in harm_fields: fields.pop(f)
		_fields = {f:fields.get(f, '.') for f in harm_fields}
		_data = self.harmonized_tuple(**_fields)
		return _data
	#-------------------------------- Harmonize ---------------------------------------
	@staticmethod
	def getSampleVAF(sample, caller, sample_ref = None):
		""" Use DP4 instead of DP
			Parameters
			----------
				sample:
				caller: {'muse', 'mutect', 'somaticsniper', 'strelka', 'varscan'}
				sample_ref: {'A', 'C', 'G', 'T'}
					Only required for strelka.
			Returns
			-------
				result: dict<>
					* 'alleles': The number of alternate alleles.
					* 'reads': The number of reads used to calculate the vaf.
					* 'vaf': The variant allele frequency.
		"""
		filter_out = False
		if caller == 'somaticsniper':
			sample_reads = sum(sample['DP4'])
			sample_alleles = sum(sample['DP4'][2:])
			sample_vaf = sample_alleles / sample_reads

		elif caller == 'mutect':
			#print(sample)
			sample_reads = sample['ALT_F1R2'] + sample['ALT_F2R1'] + sample['REF_F1R2'] + sample['REF_F2R1']
			sample_alleles = sample['ALT_F1R2'] + sample['ALT_F2R1']
			sample_vaf = sample['AF']

		elif caller == 'muse':
			sample_reads = sample['DP']
			sample_alleles = sample['AD'][1]
			sample_vaf = sample_alleles / sample_reads

		elif caller == 'strelka-snv':
			alleles = [i for i in ['A', 'C', 'G', 'T'] if i != sample_ref]
			sample_reads = sum([sample[i+'U'][1] for i in (alleles + [sample_ref])])
			sample_alleles = sum([sample[i+'U'][1] for i in alleles])
			if sample_reads == 0: sample_vaf = 0
			else:
				sample_vaf = sample_alleles / sample_reads

		elif 'varscan' in caller:
			sample_reads = sample['DP']
			sample_alleles = sample['AD']
			sample_vaf = float(sample['FREQ'].strip('%')) / 100
		else:
			sample_reads = 0
			sample_alleles = 0
			sample_vaf = 0
			print("WARNING (getSampleVAF): {0} is not a caller!".format(caller))

		result = {
			'ALLELES': sample_alleles,
			'READS': sample_reads,
			'VAF': float("{0:.5f}".format(sample_vaf)),
		}
		return result

	def _harmonizeRecord(self, record, caller):
		new_samples = list()
		for s in record.samples:
			sample_data = s.data._asdict()
			sample_vaf = self.getSampleVAF(s, caller, sample_ref = record.REF)
			sample_data = dict(list(sample_data.items()) + list(sample_vaf.items()))
			s.data = self._get_merged_data(sample_data)
			new_samples.append(s)

		new_record = vcf.model._Record(
			CHROM = record.CHROM,
			POS = record.POS,
			ID = record.ID,
			REF = record.REF,
			ALT = record.ALT,
			QUAL = record.QUAL,
			FILTER = record.FILTER,
			INFO = record.INFO,
			FORMAT = ':'.join(new_samples[0].data._fields),
			sample_indexes = record._sample_indexes,
			samples = new_samples)
		return new_record
	
	@staticmethod
	def _getHeaders(formats, infos, caller):
		new_formats = formats
		new_infos = infos
		_format_tuple = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
		#_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
		##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
		if caller == 'strelka':
			new_formats['GT'] = _format_tuple('GT', '1', 'String', 'Genotype')
		elif caller == 'varscan':
			new_formats['DP4'] = _format_tuple('DP4', '4', 'Integer', 'ref/fwd, ref/rev, var/fwd, var/rev')

		new_formats['READS'] = _format_tuple('READS', '1', 'Integer', 'Total Quality Reads')
		new_formats['VAF'] = _format_tuple('VAF', '1', 'Float', 'Variant Allele Frequency')
		new_formats['ALLELES'] = _format_tuple('ALLELES', '1', 'String', 'The number of alternate alleles')
		#new_formats['VS'] = _format_tuple('VS', '1', 'Integer', 'The Validation Status of the mutation.')

		new_infos['VS'] = _format_tuple('VS', '1', 'Integer', 'Validation Status')

		headers = {
			'formats': new_formats,
			'infos':   new_infos
		}
		return headers
	
	def _harmonizeVCFs(self, sample, variants, output_folder):
		""" Harmonizes all vcf for a given patient.
			Returns
			-------
				variants: dict<caller: path>
					A dictionary pointing to the harmonized vcf for each caller.
		"""
		for caller, filename in variants.items():
			hfile = self._harmonizeVCF(sample, caller, filename, output_folder)
			#hfile = filename
			variants[caller] = hfile
		return variants
	
	def _harmonizeVCF(self, sample, caller, filename, output_folder):
		""" Harmonizes the output VCFs from each caller.
			Parameters
			----------
				filename: string [PATH]
					The location of the output vcf. The name of the caller should be in the file's name.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.vcf
				output_folder: string [PATH]
					The folder to save the output to.
				caller: {}; defualt None
			Returns
			-------
				output_filename: string
					Path the the harmonized VCF.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Class Attributes
			----------------
				self.tag: string
					A tag to append to include when naming the file.
			
		"""        
		if caller: pass
		elif 'muse' in filename: caller = 'muse'
		elif 'mutect' in filename: caller = 'mutect'
		elif 'somaticsniper' in filename: caller = 'somaticsniper'
		elif 'strelka' in filename: caller = 'strelka'
		else: caller = 'varscan'

		source_folder, basename = os.path.split(filename)
		output_filename = "{0}_vs_{1}.{2}.harmonized.vcf".format(sample['NormalID'], sample['SampleID'], caller)
		output_filename = os.path.join(output_folder, output_filename)

		with open(filename, 'r') as input_file:
			vcf_reader = vcf.Reader(input_file)
			headers = self._getHeaders(vcf_reader.formats, vcf_reader.infos, caller)
			vcf_reader.infos   = headers['infos']
			vcf_reader.formats = headers['formats']
			#vcf_reader.formats = self._getFormats(vcf_reader.formats, caller)
			with open(output_filename, 'w') as output_file:
				vcf_writer = vcf.Writer(output_file, vcf_reader)
				for index, record in enumerate(vcf_reader):
					new_record = self._harmonizeRecord(record, caller)
					vcf_writer.write_record(new_record)
		return output_filename
	#---------------------------------- Merge -----------------------------------------
	#Validate
	def _filterOut(self, record, caller = None):
		filter_out = not record.is_snp
		if record.FILTER is None:
			return filter_out
		filters = [i for i in record.FILTER if i not in ['PASS', 'Tier1', 'Tier2', 'Tier3', 'Tier4']]
		#Caller-specific filters
		if caller == 'muse':
			if len(record.FILTER) != 0 and record.FILTER[0] not in filters:
				filter_out = True
		elif caller == 'mutect':
			if len(filters) != 0: filter_out = True
		else:
			#General caller filter
			if len(filters) != 0: filter_out = True
		if record.INFO['set'] == 'FilteredInAll': filter_out = True

		return filter_out

	def _getValidationStatus(self, record, caller = None, validation_type = 'Intersection'):
		"""
			Should use better Criteria.
			
			Use the intersection? (n = 4 or 5)
			Compare against dbSNP, COSMIC?
			Lower score if filtered by one caller but passed y others?
		"""

		normal = [i for i in record.samples if i.sample == 'NORMAL'][0]
		tumor  = [i for i in record.samples if i.sample == 'TUMOR'][0]
		if validation_type == 'VAF':
			tumor_vaf = tumor.data.VAF
			normal_vaf= normal.data.VAF

			#Filter out variants that were rejected by a caller
			#This is only needed for variants in the callsets of a single caller
			filter_out = self._filterOut(record, caller)
			
			#Validate variants according to VAF status.
			#This is only for testing purposes, a better version should be used later.
			_somatic_vaf = (tumor_vaf >= 0.08 and normal_vaf < 0.03) or (tumor_vaf < 0.08 and normal_vaf == 0.0)
			if _somatic_vaf and not filter_out:
				validation = 'Somatic'
			elif (tumor_vaf < 0.08 and (normal_vaf > 0.0 and normal_vaf < 0.03)) or filter_out:
				validation = 'Unknown'
			else: validation = 'Non-Somatic'
		elif validation_type == 'Intersection':
			callset = record.INFO
			num_callers = len(callset['set'].split('_'))
			if num_callers >= 4:
				validation = 'Somatic'
			else:
				validation = 'Non-Somatic'
			#print(num_callers, callset['set'])
				
		if validation == 'Non-Somatic': validation = 0
		elif validation == 'Somatic': validation = 1
		else: validation = 2
		
		position = "{0}:{1}".format(record.CHROM, record.POS)
		validation = self.variant_validation_status.get(position, 0)

		return validation
	
	#Merge
	def merge_vcfs(self,sample, options, variants, output_folder):
		""" Uses GATK CombineVariants to merge the calls from each caller into a single file.
			Parameters
			----------
				variants: dict<caller, path>
					A dictionary linkng each caller to its harmonized output.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Returns
			-------
				Output_filename: string
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.merged.vcf
		"""
		if os.name == 'nt':
			gatk = "C:\\Users\\Deitrickc\\Downloads\\Genomic Programs\\GenomeAnalysisTK-3.7\\GenomeAnalysisTK.jar"
			reference = "G:\\Pipeline Files\\Reference\\GRCh38.d1.vd1.fa"
		else:
			gatk = options['Programs']['GATK']
			reference = options['Reference Files']['reference genome']
		
		#Check GATK
		if not os.path.isfile(gatk):
			print("ERROR: Could not locate GATK at ", gatk)
		if not os.path.isfile(reference):
			print("ERROR: Could not locate reference at ", reference)

		output_file = "{0}_vs_{1}.merged.vcf".format(sample['NormalID'], sample['SampleID'])
		output_file = os.path.join(output_folder, output_file)
		
		order = "mutect,varscan,strelka,muse,somaticsniper" #ordered by VAF confidence
		#order = ",".join([i for i in order.split(',') if i in variants])
		
		template = "--variant:{0} {1}"
		templates = ' \\\n'.join([template.format(c, f) for c, f in variants.items()])
		
		command = """java -jar "{gatk}" \
			-T CombineVariants \
			-R "{reference}" \
			--variant:muse "{muse}" \
			--variant:mutect "{mutect}" \
			--variant:varscan "{varscan}" \
			--variant:somaticsniper "{ss}" \
			--variant:strelka "{strelka}" \
			-o "{output}" \
			-genotypeMergeOptions PRIORITIZE \
			-priority {rod}"""
		command = command.format(
				gatk = gatk,
				reference = reference,
				muse = variants['muse'],
				mutect = variants['mutect'],
				ss = variants['somaticsniper'],
				strelka = variants['strelka-snv'],
				varscan = variants['varscan-snv'],
				variants = templates,
				rod = order,
				output = output_file)

		if not os.path.isfile(output_file):
			print(command)
			os.system(command)
		return output_file
		
	def _getValidatedRecord(self, record):
		#Modify the INFO field to add the validation status and format the 'sets' field
		#infotuple = collections.namedtuple('INFO', sorted(record.INFO) + ['VS'])
		
		
		new_info = record.INFO
		if new_info['set'] == 'Intersection':
			caller_sets = sorted(['muse', 'mutect', 'somaticsniper', 'strelka', 'varscan'])
		else:
			caller_sets = sorted(i for i in new_info['set'].split('-') if 'filter' not in i)
		caller_sets = '_'.join(caller_sets)
		new_info['set'] = caller_sets
		new_info['VS'] = self._getValidationStatus(record)
		#new_info = infotuple(**new_info)

		new_record = vcf.model._Record(
			CHROM = record.CHROM,
			POS = record.POS,
			ID = record.ID,
			REF = record.REF,
			ALT = record.ALT,
			QUAL = record.QUAL,
			FILTER = record.FILTER,
			INFO = new_info,#record.INFO,
			FORMAT = record.FORMAT,
			sample_indexes = record._sample_indexes,
			samples = record.samples)
		
		return new_record
	 
	def validate_variants(self, sample, merged_vcf, output_folder):
		""" Adds a field to the merged vcf file indicating whether a given mutation is present in the truth set.

			Returns
			-------
				output_filename: string [PATH]
				Format: {normalID}_vs_{tumorID}.merged.validated.vcf
		"""
		output_filename = "{0}_vs_{1}.merged.validated.vcf".format(sample['NormalID'], sample['SampleID'])
		output_filename = os.path.join(output_folder, output_filename)
		#source_folder, basename = os.path.split(merged_vcf)
		#basename = basename.split('.')[0] + self.tag + '.validated.vcf'
		#output_filename = os.path.join(source_folder, basename)
		
		with open(merged_vcf, 'r') as file1:
			vcf_reader = vcf.Reader(file1)
			with open(output_filename, "w") as output:
				vcf_writer = vcf.Writer(output, vcf_reader)
				
				for record in vcf_reader:
					#validation_status = self._getValidationStatus(record)
					new_record = self._getValidatedRecord(record)
					vcf_writer.write_record(new_record)
		return output_filename
			
	def toTable(self, merged_vcf, output_folder):
		print("Merged VCF FIlename: ", merged_vcf)
		table = list()
		stats = collections.defaultdict(list)
		with open(merged_vcf, 'r') as file1:
			vcf_reader = vcf.Reader(file1)
			for record in vcf_reader:
				normal = [s for s in record.samples if s.sample == 'NORMAL'][0].data._asdict()
				tumor  = [s for s in record.samples if s.sample == 'TUMOR'][0].data._asdict()
				variant_allele = record.ALT[0]
				if variant_allele is None or len(variant_allele) > 1: continue
				row = {
					"Patient": "TCGA-2H-A9GF",
					"Validation_Status": record.INFO['VS'],
					"Mutation_Set": record.INFO['set'],
					"Reference_Allele": record.REF,
					"Variant_Allele": variant_allele,
					"Depth_Tumor": tumor['READS'],
					"Depth_Normal": normal['READS'],
					"Vaf_Tumor": tumor['VAF']*100,
					"Vaf_Normal": normal['VAF']*100
				}
				
				if row['Mutation_Set'] == 'FilteredInAll': continue
				elif row['Validation_Status'] == 2: continue
				
				stats[row['Mutation_Set']].append(row['Validation_Status'])
				table.append(row)
		print("Validation per set:")

		#for key, series in sorted(stats.items()):
			
		#	_total_validated = series.count(1)
		#	_total_detected = len(series)
		#	_ratio = _total_validated / _total_detected
		#	print("{0:<45}\t{1}\t{2:.1%}".format(key, _total_detected, _ratio))

		basename = os.path.splitext(os.path.basename(merged_vcf))[0] + '.table.tsv'
		table_filename = os.path.join(output_folder, basename)
		header = sorted(table[0].keys())
		with open(table_filename, 'w') as outputfile:
			writer = csv.DictWriter(outputfile, delimiter = '\t', fieldnames = header)
			writer.writeheader()
			writer.writerows(table)
		return stats
	@staticmethod
	def _count_validated_variants(stats, filename):
		"""
			Parameters
			----------
				stats: dict<list<>>
					A dictionary mapping callset names to a list of callset varaints validation data.
					Ex. MuSE: [0,1,1,1,0,0,0,0,1,00,1,0,1]
		"""

		with open(filename, 'w') as file1:
			file1.write("Callset\ttotal calls\ttotal validated\tpercent validated\n")
			for key, series in stats.items():
				_total_validated = series.count(1)
				_total_detected = len(series)
				_ratio = _total_validated / _total_detected
				file1.write("{0:<45}\t{1}\t{2}\t{3:.1%}\n".format(key, _total_detected, _total_validated, _ratio))
	def expandTable(table_file, callers, caller_sets):
		""" Adds and initializes columns to a dataframe based on the passed variable list.
		Parameters
		----------
			df: pandas.DataFrame
				The dataframe to modify
			ggname: string
				The name of an existing column in the df.
			filters: list<str>
		
		"""
		#The subtype is just the reference and alternate alleles concatenated
		for line in table:
			row = line
			row['subtype'] = row['Reference_Allele'] + row['Variant_Allele']
			for caller in callers:
				row['caller'] = int(caller in row['Mutation_Set'])
			for caller_set in caller_sets:
				row[caller_set] = int(row['Mutation_Set'] == caller_set)

		return df

class Truthset:
	"""
		Input
		-----
			Intersection


		Output
		------
			The output is a vcf file of all sites that pass any truthset filters. the VAF for each
			site is included.

			Folder: PIPELINE_FOLDER/PROCESSED_VCFS_FOLDER/truthset/
			Filename: [NORMAL]_vs_[TUMOR].[Training type].['snv'/'indel'].truthset.vcf (includes index file)
		Description
		-----------

	"""
	def __init__(self, samples, options, training_type, **kwargs):
		""" Generates a truthset based on one or more samples.
			Option 1: intersection of all 5 callers.
			Option 2: Dream-SEQ data
			Option 3: RNA-seq?
			Parameters
			----------
				samples: dict<>
					A dictionary containing information for the patient.
				options: dict<>
				variants: dict<caller, filename>
				training_type: {'intersection'}

				**kwargs:
					* 'intersection': int; default 5
						number of callers to pass validation when 'training_type' is 'intersection'
		"""
		kwargs['intersection'] = kwargs.get('intersection', 5)
		self.training_type = training_type
		outputs = list()
		for sample in samples:
			sample_truthset = self._per_sample(sample, options, training_type, **kwargs)
			outputs.append(sample_truthset)

		if len(outputs) == 1:
			sample_truthset = outputs[0]
		else:
			sample_truthset = self._combineTruthsets(outputs, options)
		#self.truthset = sample_truthset['truthset']
		self.filename = sample_truthset['filename']
		self.indel_filename = sample_truthset['filename-indel']
		self.snv_filename = sample_truthset['filename-snv']

		#pprint(sample_truthset)


	def __init__OBS(self, sample, options, training_type, **kwargs):
		""" Generates a truthset based on the passed variants/
			Option 1: intersection of all 5 callers.
			Option 2: Dream-SEQ data
			Option 3: RNA-seq?
			Parameters
			----------
				sample: dict<>
				options: dict<>
				variants: dict<caller, filename>
				training_type: {'intersection'}

				**kwargs:
					* 'n': int; default 5
						number of callers to pass validation when 'training_type' is 'intersection'
		"""
		#Filename format: TCGA-2H-A9GF-11A_vs_TCGA-2H-A9GF-01A.truthset.RNA-seq.txt
		
		#Define the output file
		print("Truthset(training_type = {0})".format(training_type))
		self.training_type = training_type
		input_vcf, output_vcf, self.snv_filename, self.indel_filename = self._get_vcf_files(
			sample, options, training_type)
		self.filename = output_vcf
		#Select truthset source
		if not os.path.exists(output_vcf):
			self.truthset = self._generate_truthset(input_vcf, output_vcf, training_type)
		else:
			self.truthset = vcf.Reader(open(output_vcf, 'r'))

		self._split_vcf()
	
	def __call__(self, sample, chrom, pos):
		"""
			Parameters
			----------
				sample: tumor barcode
				chrom: int
				pos: int
		"""
		if self.isTable:
			rows = self.truthset.extract()
		else:
			rows = [i for i in self.truthset if i['sample'] == sample and i['chrom'] == chrom and i['position'] == pos]

			response = {row['validation method']: row['validation status'] for row in rows}
		
		for i in ['RNA-seq', 'Intersection', 'VAF']:
			if i not in response.keys(): response[i] = 0

		return response

	def _per_sample(self, sample, options, training_type, **kwargs):
		""" Generates individual truthsets per sample.
			Available Keyword Arguments
			---------------------------
				'intersection': number of callers to count as the intersection.
				'variant_type': If not 'snv' or 'indel', the file will be split into
					an 'snv' and 'indel' file.
		"""
		vcf_filenames = self._get_vcf_files(sample, options, training_type)
		output_vcf = vcf_filenames['output']
		snv_filename = vcf_filenames['snv']
		indel_filename = vcf_filenames['indel']
		#self.filename = output_vcf
		#Select truthset source
		output_folder = options['Pipeline Options']['processed vcf folder']
		output_folder = os.path.join(output_folder, 'merged_vcfs', sample['PatientID'])
		checkdir(output_folder)
		
		if training_type == 'RNA':
			variants = GetVariantList(sample, options, 'RNA-seq')

		else:
			variants = GetVariantList(sample, options, 'DNA-seq')
		combined_variants = CombineVariants(sample, options, output_folder, variants = variants)
		
		#combined_variants = CombineVariants(sample, options, output_folder)
		snv_variants = combined_variants.snvs
		indel_variants = combined_variants.indels
		
		print("Output VCF: ", output_vcf)
		#if not os.path.exists(output_vcf) or True:
		snv_truthset = self._generate_truthset(snv_variants, snv_filename, training_type, **kwargs)
		indel_truthset = self._generate_truthset(indel_variants, indel_filename, training_type, **kwargs)
		#else:
		#	truthset = vcf.Reader(open(output_vcf, 'r'))

		#if kwargs['variant_type'] == 'all':
		#	self._split_vcf(output_vcf, snv_filename, indel_filename)

		result = {
			'PatientID': sample['PatientID'],
			#'truthset': truthset,
			'filename': output_vcf,
			'filename-indel': indel_filename,
			'filename-snv': snv_filename
		}
		return result

	def _split_vcf(self, filename, snv_filename, indel_filename):
		""" Splits the full vcf file into snvs and indels """
		reader = vcf.Reader(open(filename, 'r'))
		snv_writer = vcf.Writer(open(snv_filename, 'w'), reader)
		indel_writer = vcf.Writer(open(indel_filename, 'w'), reader)

		for record in reader:
			if record.is_snp:
				snv_writer.write_record(record)
			elif record.is_indel:
				indel_writer.write_record(record)

	@staticmethod
	def _get_vcf_files(sample, options, training_type):
		""" Keyword Arguments
			-----------------
				variant_type
		"""
		output_vcf = os.path.join(
			options['Pipeline Options']['processed vcf folder'],
			"truthset",
			"{0}_vs_{1}.{2}.truthset.vcf".format(sample['NormalID'], sample['SampleID'], training_type))

		snv_vcf = os.path.join(
			options['Pipeline Options']['processed vcf folder'],
			"truthset",
			"{0}_vs_{1}.{2}.snv.truthset.vcf".format(sample['NormalID'], sample['SampleID'], training_type))
		
		indel_vcf = os.path.join(
			options['Pipeline Options']['processed vcf folder'],
			"truthset",
			"{0}_vs_{1}.{2}.indel.truthset.vcf".format(sample['NormalID'], sample['SampleID'], training_type))

		checkdir(os.path.dirname(output_vcf))
		"""
		if training_type in {'Intersection', 'VAF'}:
			#Both types require the output from GATK Combine Variants.
			input_vcf = os.path.join(
				options['Pipeline Options']['processed vcf folder'], 
				"merged_vcfs", 
				"{0}_vs_{1}.merged.vcf".format(sample['NormalID'], sample['SampleID']))
		elif training_type == 'RNA-seq':
			input_vcf = os.path.join(
				options['Pipeline Options']['somatic pipeline folder'],
				sample['PatientID'],
				'HaplotypeCaller', 'RNA-seq',
				'TCGA-2H-A9GF-01A.RNA.raw_snps_indels.vcf')
		"""
		result = {
			'output': output_vcf,
			'snv': snv_vcf,
			'indel': indel_vcf
		}
		return result
	@staticmethod
	def _getSampleVAF(record, caller = None):
		""" Use DP4 instead of DP
			Parameters
			----------
				record: from pyVCF reader
			Returns
			-------
				result: dict<>
				* {NORMAL, SAMPLE}:
					* 'alleles': THe number of alternate alleles.
					* 'reads': The number of reads used to calculate the vaf.
					* 'vaf': The variant allele frequency.
		"""
		response = dict()
		#record = record.samples
		for sample in record.samples:
			sample_fields = sample.data._asdict()
			if caller is None:
				if 'DP4' in sample_fields: 			caller = 'somaticsniper'
				elif 'ALT_F1R2' in sample_fields: 	caller = 'mutect'
				elif 'FREQ' in sample_fields: 		caller = 'varscan'
				elif 'DP' in sample_fields: 		caller = 'muse'
				else: 								caller = 'strelka'
			else:
				caller = caller.lower()
			#print(caller)
			#pprint(sample_fields)
			if caller == 'somaticsniper':
				
				sample_reads = sum(sample_fields['DP4'])
				sample_alleles = sum(sample_fields['DP4'][2:])
				sample_vaf = sample_alleles / sample_reads

			elif caller == 'mutect':
				#print(sample)
				sample_reads = sample_fields['ALT_F1R2'] + sample_fields['ALT_F2R1'] + sample_fields['REF_F1R2'] + sample_fields['REF_F2R1']
				sample_alleles = sample_fields['ALT_F1R2'] + sample_fields['ALT_F2R1']
				sample_vaf = sample_fields['AF']

			elif caller == 'muse':
				sample_reads = sample_fields['DP']
				sample_alleles = sample_fields['AD'][1]
				sample_vaf = sample_alleles / sample_reads

			elif caller == 'strelka':
				sample_ref = record.REF
				alleles = [i for i in ['A', 'C', 'G', 'T'] if i != sample_ref]
				sample_reads = sum([sample_fields[i+'U'][1] for i in (alleles + [sample_ref])])
				sample_alleles = sum([sample_fields[i+'U'][1] for i in alleles])
				if sample_reads == 0: sample_vaf = 0
				else:
					sample_vaf = sample_alleles / sample_reads

			elif caller == 'varscan':
				sample_reads = sample_fields['DP']
				sample_alleles = sample_fields['AD']
				sample_vaf = float(sample_fields['FREQ'].strip('%')) / 100
			else:
				print("WARNING (getSampleVAF): {0} is not a caller!".format(caller))

			result = {
				'reads': sample_reads,
				'alleles': sample_alleles,
				'vaf': sample_vaf
			}
			response[sample.sample] = result

		return response

	def _generate_truthset(self, input_vcf, output_vcf, training_type, **kwargs):
		""" Generates a truthset
			Parameters
			----------
				sample
				options
				training_type: {'RNA-seq', 'VAF', 'Intersection'}
		 """
		with open(input_vcf, 'r') as input_file:
			reader = vcf.Reader(input_file)
			with open(output_vcf, 'w') as output_file:
				writer = vcf.Writer(output_file, reader)

				for index, record in enumerate(reader):
					recordStatus = self._getRecordValidationStatus(record, training_type, **kwargs)
					isValid = recordStatus['validation status']
					
					if isValid:
						writer.write_record(record)
					#lif record.is_indel and isFirstIndel:
					#	writer.write_record(record)
					#	isFirstIndel = False
		#reader = vcf.Reader(open(output_vcf, 'r'))
		reader = None
		return reader

	def _getRecordValidationStatus(self, record, training_type, **kwargs):
		""" Determines if a given record is within the truthset.
			Returns
			-------
				dict<>
					* 'chrom': 
					* 'position': 
					* 'validation method': 
					* 'validation status': 
		"""
		if training_type == 'Intersection':
			recordStatus = self._from_intersection(record, **kwargs)
		elif training_type == 'RNA-seq':
			#Assume the record is from the vcf from the RNA-seq pipeline.
			#Assume all RNA-seq variants are true variants.
			recordStatus = self._from_RNA(record)
		elif training_type == 'VAF':
			recordStatus = self._from_VAF(record)

		return recordStatus

	def _combineTruthsets(self, samples, options):
		""" Combines the truthsets of several samples into a single file.
			Parameters
			----------
				samples: list<dict<>>
					A list of outputs from self._per_sample()
						* 'filename-indel'
						* 'filename-snv'
		"""
		gatk_program = options['Programs']['GATK']
		program = options['Programs']['Picard']
		reference = options['Reference Files']['reference genome']

		_folderName = lambda s: os.path.join(options['Pipeline Options']['processed vcf folder'], 'truthset', s)
		sampleIds = sorted(i['PatientID'] for i in samples)
		sampleIds = [i.split('-')[-1] for i in sampleIds] #Extracts patient id from the barcode. ex. TCGA-2H-A9GF -> A9GF
		sampleIds = ",".join(sampleIds)
		
		filename 		= _folderName("{0}.{1}.merged_truthset.vcf".format(sampleIds, self.training_type))
		indel_filename 	= _folderName("{0}.{1}.merged_truthset.indel.vcf".format(sampleIds, self.training_type))
		snv_filename 	= _folderName("{0}.{1}.merged_truthset.snv.vcf".format(sampleIds, self.training_type))

		#variants       = "--variant " + " --variant ".join([i['filename'] for i in samples])
		snv_variants   = "--variant " + " --variant ".join(i['filename-snv'] for i in samples) + ' \\'
		indel_variants = "--variant " + " --variant ".join([i['filename-indel'] for i in samples]) + '\\'
		snv_variants   = " ".join(["I={0}".format(i['filename-snv']) for i in samples])
		indel_variants   = " ".join(["I={0}".format(i['filename-indel']) for i in samples])
		gatk_base_command = """ java -cp {program} org.broadinstitute.gatk.tools.CatVariants \
				    -R {reference} \
				    {variants}
				    -out {output}"""
		picard_base_command = """java -jar {program} SortVcf \
			{variants} \
			O={output}"""
		snv_command = picard_base_command.format(
			program = program,
			reference = reference,
			variants = snv_variants,
			output = snv_filename)
		"""
		java -jar picard.jar SortVcf \
      I=vcf_1.vcf \
      I=vcf_2.vcf \
      O=sorted.vcf
		"""
		indel_command = picard_base_command.format(
			program = program,
			reference = reference,
			variants = indel_variants,
			output = indel_filename)

		print(snv_command)
		os.system(snv_command)
		os.system(indel_command)
		result = {
			'PatientID': sampleIds,
			'filename': filename,
			'filename-indel': indel_filename,
			'filename-snv': snv_filename
		}

		return result

	def _combineVariants(self, sample, options):
		""" Generates a merged vcf file, if one is not already present."""
	def _generate_truthsetOBS(self, sample, options, training_type):
		print("Truthset._generate_table({0}, {1}, {2})".format(type(sample), type(options), training_type))
		table = list()
		training_files = dict()
		if training_type == 'RNA-seq' or training_type is None:
			filename = os.path.join(
				options['Pipeline Options']['somatic pipeline folder'],
				sample['PatientID'],
				'HaplotypeCaller', 'RNA-seq',
				'TCGA-2H-A9GF-01A.RNA.raw_snps_indels.vcf')
			training_files['RNA-seq'] = [filename]
		
		if training_type == 'Intersection':
			filename = os.path.join(
				options['Pipeline Options']['processed vcf folder'], 
				"merged_vcfs", 
				"{0}_vs_{1}.merged.vcf".format(sample['NormalID'], sample['SampleID']))
			training_files['Intersection'] = [filename]
		
		if training_type == 'VAF' or training_type is None:
			training_files['VAF'] = GetVariantList(sample, options)
		
		pprint(training_files)
		for t_type, filenames in training_files.items():
			for filename in filenames:
				print(t_type, filename)
				if isinstance(filenames, dict):
					caller, filename = filename, filenames[filename]
				else: caller = None
				with open(filename, 'r') as vcf_file:
					vcf_reader = vcf.Reader(fsock = vcf_file)
					for index, record in enumerate(vcf_reader):
						#print(record)
						if t_type == 'Intersection':
							row = self._from_intersection(record)
						elif t_type == 'RNA-seq':
							row = self._from_RNA(record)
						elif t_type == 'VAF':
							row = self._from_VAF(record, caller = caller)
						else:
							print(t_type, " is not a valid training type!")
						row['sample'] = sample['SampleID']
						table.append(row)
		return table
	
	def _from_intersection(self, record, **kwargs):
		if '-' in record.INFO['set']:
			_separator = '-'
		else: _separator = '_'

		#GATK includes callers that filtered the variant, designated by a "filterIn{caller}" label.
		callset = [i for i in record.INFO['set'].split(_separator) if 'filter' not in i.lower()]

		#MuSE and Somaticsniper are snp-only, and should not be counted in the intersection of indel callsets
		if record.is_indel:
			num_callers_in_intersection = 3
			#num_callers_in_intersection = kwargs['indel_intersection']
		else:
			num_callers_in_intersection = 5
			#num_callers_in_intersection = kwargs['snv_intersection']


		_is_intersection = len(callset) == 1 and callset[0] == 'Intersection'
		_is_in_n_sets = int(len(callset) >= num_callers_in_intersection)
		validation_status = int(_is_intersection or _is_in_n_sets)
		
		#if _is_intersection: validation_status = 5
		#else: validation_status = len(callset)


		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'Intersection',
			'validation status': validation_status
		}
		return row

	def _from_dream(self, sample, filename):
		pass
	@staticmethod
	def _from_RNA(record):
		filterOut = False
		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'RNA-seq',
			'validation status': int(filterOut == False)
		}
		return row

	def _from_VAF(self, record, caller):
		sample_vaf = self._getSampleVAF(record, caller)
		normal_vaf = sample_vaf['NORMAL']['vaf']
		tumor_vaf = sample_vaf['TUMOR']['vaf']


		#Filter out variants that were rejected by a caller
		#This is only needed for variants in the callsets of a single caller
		#filter_out = self._filterOut(record, caller)
		filter_out = False
		#Validate variants according to VAF status.
		#This is only for testing purposes, a better version should be used later.
		_somatic_vaf = (tumor_vaf >= 0.08 and normal_vaf < 0.03) or (tumor_vaf < 0.08 and normal_vaf == 0.0)
		if _somatic_vaf and not filter_out:
			validation_status = 1 #Somatic
		elif (tumor_vaf < 0.08 and (normal_vaf > 0.0 and normal_vaf < 0.03)) or filter_out:
			validation_status = 0 #Non-Somatic
		else: validation_status = 2 #'Unknown'

		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'VAF',
			'validation status': validation_status
		}
		return row
	
	def to_vcf(table, filename):
		pass
	@staticmethod
	def save(table, filename):
		table = sorted(table, key = lambda s: (s['sample'], s['chrom'], s['position']))
		with open(filename, 'w', newline = '') as file1:
			#reader = csv.DictReader(file1, delimiter = '\t')
			writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = ['sample', 'chrom', 'position', 'validation method', 'validation status'])
			writer.writeheader()
			writer.writerows(table)
		return filename

def GetVariantList(sample, options = None, folder = None, seq_type = 'DNA-seq'):
	""" Returns a map of each caller's indel and snv outputs.
		Parameters
		----------
			folder: string
				The folder containing the somatic variants for all samples.
			barcode: string
				The TCGA barcode for the case
	"""
	#spf = options['Pipeline Options']['somatic pipeline folder']
	if folder is None:
		spf = SOMATIC_PIPELINE_FOLDER
	else:
		spf = folder
	vcf_folder = os.path.join(spf, barcode)
	if seq_type == 'DNA-seq':


		patientID = sample['PatientID']
		normalID  = sample['NormalID']
		tumorID   = sample['SampleID']

		variants   = {
			#os.path.join(options['output']['Bambino'].format(patient = patientID),"{0}_vs_{1}.bambino.vcf".format(normalID, tumorID)),
			#os.path.join(options['output']['Haplotypecaller'].format(patient = patientID), "{0}_vs_{1}.raw.snps.indels.vcf".join(normalID, tumorID)),
			'muse':          os.path.join(vcf_folder, "MuSE", 				"{0}_vs_{1}.Muse.vcf".format(normalID, tumorID)),
			'mutect':        os.path.join(vcf_folder, "MuTect2", 			"{0}_vs_{1}.mutect2.vcf".format(normalID, tumorID)),
			'somaticsniper': os.path.join(vcf_folder, "SomaticSniper", 		"{0}_vs_{1}.somaticsniper.hq.vcf".format(normalID, tumorID)),
			'strelka-indel': os.path.join(vcf_folder, "Strelka", "results", "{0}_vs_{1}.passed.somatic.indels.vcf.strelka.vcf".format(normalID, tumorID)),
			'strelka-snv':   os.path.join(vcf_folder, "Strelka", "results", "{0}_vs_{1}.passed.somatic.snvs.vcf.strelka.vcf".format(normalID, tumorID)),
			'varscan-indel': os.path.join(vcf_folder, "Varscan", 			"{0}_vs_{1}.raw.indel.vcf".format(normalID, tumorID)),
			'varscan-snv':   os.path.join(vcf_folder, "Varscan", 			"{0}_vs_{1}.raw.snp.Somatic.hc.vcf".format(normalID, tumorID))
		}
	else:

		variants = {
			'haplotypecaller': os.path.join(vcf_folder, "HaplotypeCaller", "RNA-seq", "{0}.RNA.raw_snps_indels.vcf".format(sample['SampleID']))
		}

	return variants

def VariantEffectPredictor(sample, options, merged_variants):
	""" VariantEffectPredictor uses terminology as defined by the Sequence Ontology Project
	"""
	#perl variant_effect_predictor.pl --cache -i input.txt -o output.txt	
	print("Running VEP...", flush = True)
	pprint(merged_variants)
	
	output_folder = os.path.join(options['Pipeline Options']['processed vcf folder'], "vep_vcfs")
	reference = options['Reference Files']['reference genome']
	program = options['Programs']['varianteffectpredictor']
	
	checkdir(output_folder)

	annotated_variants = dict()
	for caller, source in merged_variants.items():
		_, fn = os.path.split(source)
		destination = os.path.join(output_folder, os.path.splitext(fn)[0] + '.vep.vcf')

		command = """perl {vep} \
			--input_file {inputfile} \
			--output_file {output} \
			--fasta {reference} \
			--species homo_sapiens \
			--assembly GRCh38 \
			--format vcf \
			--cache \
			--html \
			--symbol \
			--biotype \
			--total_length \
			--numbers \
			--fields Consequence,Codons,Aminoacids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
			--vcf""".format(
				vep = program,
				inputfile = source,
				output = destination,
				reference = reference)
		#print(command)
		if not os.path.exists(destination):
			Terminal(command, show_output = True)

		annotated_variants[caller] = destination

	return annotated_variants

class toMAF:
	def __init__(self, sample, options, input_vcf, truthset = None):
		""" Converts a VCF file to and annotated MAF file. Include additional fields such as VAF.
			Additional Columns
			------------------
				VAF: The variant allele frequency

		"""
		if truthset is None:
			truthset = Truthset(sample, options)

		maf_file = self.vcftomaf(sample, options, input_vcf)
		maf_file = self.modifyMAF(sample, options, maf_file, truthset)


		output_folder = options['Pipeline Options']['maf folder']
		output_filename = os.path.join(
			output_folder, 
			sample['PatientID'],
			os.path.splitext(os.path.basename(input_vcf))[0] + '.raw.maf')

		checkdir(os.path.dirname(output_filename))
		print("Saving to ", output_filename)
		
		with open(output_filename, 'w', newline = '') as file1:
			file1.write("#version 2.4\n")
			fieldnames = "Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome	Start_Position	End_Position	Strand	Variant_Classification	Variant_Type	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	dbSNP_RS	dbSNP_Val_Status	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	Match_Norm_Seq_Allele1	Match_Norm_Seq_Allele2	Tumor_Validation_Allele1	Tumor_Validation_Allele2	Match_Norm_Validation_Allele1	Match_Norm_Validation_Allele2	Verification_Status	Validation_Status	Mutation_Status	Sequencing_Phase	Sequence_Source	Validation_Method	Score	BAM_File	Sequencer	Tumor_Sample_UUID	Matched_Norm_Sample_UUID	HGVSc	HGVSp	HGVSp_Short	Transcript_ID	Exon_Number	t_depth	t_ref_count	t_alt_count	n_depth	n_ref_count	n_alt_count	all_effects	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	ALLELE_NUM	DISTANCE	STRAND_VEP	SYMBOL	SYMBOL_SOURCE	HGNC_ID	BIOTYPE	CANONICAL	CCDS	ENSP	SWISSPROT	TREMBL	UNIPARC	RefSeq	SIFT	PolyPhen	EXON	INTRON	DOMAINS	GMAF	AFR_MAF	AMR_MAF	ASN_MAF	EAS_MAF	EUR_MAF	SAS_MAF	AA_MAF	EA_MAF	CLIN_SIG	SOMATIC	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	IMPACT	PICK	VARIANT_CLASS	TSL	HGVS_OFFSET	PHENO	MINIMISED	ExAC_AF	ExAC_AF_AFR	ExAC_AF_AMR	ExAC_AF_EAS	ExAC_AF_FIN	ExAC_AF_NFE	ExAC_AF_OTH	ExAC_AF_SAS	GENE_PHENO	FILTER	flanking_bps	variant_id	variant_qual	ExAC_AF_Adj	ExAC_AC_AN_Adj	ExAC_AC_AN	ExAC_AC_AN_AFR	ExAC_AC_AN_AMR	ExAC_AC_AN_EAS	ExAC_AC_AN_FIN	ExAC_AC_AN_NFE	ExAC_AC_AN_OTH	ExAC_AC_AN_SAS	ExAC_FILTER".split('\t')
			fieldnames += sorted(i for i in maf_file[0].keys() if i not in fieldnames)
			print(fieldnames)
			writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = fieldnames)
			writer.writeheader()
			writer.writerows(maf_file)
		self.filename = output_filename

	@staticmethod
	def modifyMAF(sample, options, input_maf, truthset):
		""" Adds relevant information to the maf files
			Additional Columns
			------------------
				Tumor_Sample_Barcode: SampleID
				Matched_Norm_Sample_Barcode: NormalID
				tumor_bam_uuid: SampleUUID
				normal_bam_uuid: NormalUUID

				Validation_Status_VAF: Validation status determined from VAF
					0: Non-Somatic
					1: Somatic
					2: Unknown
				Validation_Status_RNA-seq:
					0: Not present in RNA-seq
					1: Present in RNA-seq
				Validation_Status_Intersect:
					0-5: Present in n callers
				Filter_Status: Whether the variant passes all required filters
				Filter_Score: float
					Probability of a variant being somatic
				Callset:
				
			Modifies Columns
			----------------
				dbSNP_Val_Status: Whether the variant is in dbSNP
				Sequencer:
		"""
		print("toMAF.modifyMAF(sample, options, {0}, truthset)".format(input_maf))
		output_maf = list()

		with open(input_maf, 'r') as maf_file:
			maf_file.seek(0)
			next(maf_file)
			reader = csv.DictReader(maf_file, delimiter = '\t')

			for index, row in enumerate(reader):
				if index % 500 == 0: print(index)

				row['Tumor_Sample_Barcode'] 		= sample['SampleID']
				row['Matched_Norm_Sample_Barcode'] 	= sample['NormalID']
				row['tumor_bam_uuid'] 				= sample['SampleUUID']
				row['normal_bam_uuid'] 				= sample['NormalUUID']
				#row['Callset'] 						= 'MuSE'
				#pprint(row)
				validation_status = truthset(sample = sample['SampleID'], chrom = row['Chromosome'], pos = row['Start_Position'])			

				row['Validation_Status_VAF'] 		= validation_status.get('VAF')
				row['Validation_Status_RNA'] 		= validation_status.get('RNA-seq')
				row['Validation_Status_Intersect'] 	= validation_status.get('Intersection')

				row['Filter_Status'] = 0

				output_maf.append(row)

		return output_maf

	@staticmethod
	def vcftomaf(sample, options, input_vcf):
		"""
		 --input-vcf      Path to input file in VCF format
		 --output-maf     Path to output MAF file [Default: STDOUT]
		 --tmp-dir        Folder to retain intermediate VCFs after runtime [Default: Folder containing input VCF]
		 --tumor-id       Tumor_Sample_Barcode to report in the MAF [TUMOR]
		 --normal-id      Matched_Norm_Sample_Barcode to report in the MAF [NORMAL]
		 --vcf-tumor-id   Tumor sample ID used in VCF's genotype columns [--tumor-id]
		 --vcf-normal-id  Matched normal ID used in VCF's genotype columns [--normal-id]
		 --custom-enst    List of custom ENST IDs that override canonical selection
		 --vep-path       Folder containing variant_effect_predictor.pl [~/vep]
		 --vep-data       VEP's base cache/plugin directory [~/.vep]
		 --vep-forks      Number of forked processes to use when running VEP [4]
		 --buffer-size    Number of variants VEP loads at a time; Reduce this for low memory systems [5000]
		 --any-allele     When reporting co-located variants, allow mismatched variant alleles too
		 --ref-fasta      Reference FASTA file [~/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz]
		 --filter-vcf     The non-TCGA VCF from exac.broadinstitute.org [~/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz]
		 --max-filter-ac  Use tag common_variant if the filter-vcf reports a subpopulation AC higher than this [10]
		 --species        Ensembl-friendly name of species (e.g. mus_musculus for mouse) [homo_sapiens]
		 --ncbi-build     NCBI reference assembly of variants MAF (e.g. GRCm38 for mouse) [GRCh37]
		 --cache-version  Version of offline cache to use with VEP (e.g. 75, 82, 86) [Default: Installed version]
		 --maf-center     Variant calling center to report in MAF [.]
		 --retain-info    Comma-delimited names of INFO fields to retain as extra columns in MAF []
		 --min-hom-vaf    If GT undefined in VCF, minimum allele fraction to call a variant homozygous [0.7]
		 --remap-chain    Chain file to remap variants to a different assembly before running VEP
		 --help           Print a brief help message and quit
		 --man            Print the detailed manual

		"""
		print("Converting VCF to MAF...", flush = True)
		vcftomaf_script = options['Programs']['vcf2maf']
		reference = options['Reference Files']['reference genome']
		vep_location = os.path.dirname(options['Programs']['varianteffectpredictor'])
		#input_vcf = merged_variants
		
		output_folder = options['Pipeline Options']['maf folder']
		output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_vcf))[0] + ".maf")

		if isinstance(input_vcf, dict): input_vcf = list(input_vcf.values())[0]
		#perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --tumor-id WD1309 --normal-id NB1308

		command = """perl {script} \
			--input-vcf {vcf} \
			--output-maf {maf} \
			--tumor-id {tumor} \
			--normal-id {normal} \
			--vcf-tumor-id TUMOR \
			--vcf-normal-id NORMAL \
			--buffer-size 100 \
			--ref-fasta {reference} \
			--ncbi-build {ncbi} \
			--vep-path {vep}""".format(
				script = vcftomaf_script,
				vcf = input_vcf,
				maf = output_file,
				tumor = sample['SampleID'],
				normal = sample['NormalID'],
				reference = reference,
				ncbi = 'GRCh38',
				vep = vep_location)
		#print(command)
		Terminal(command)

		return output_file

	def filterMAF(sample, options, input_maf, caller):
		"""
			Filter Criteria
			---------------
				dbSNP: filter out variants present in dbSNP
				Indel Filtering: filter out variants in known indels
				Germline/normal filtering: filter out variants present in a panel of normals. (derive from GDC MAFs?)
				PON Filtering
				
		"""
		return row

def SomaticSeq(sample_list, options, training_type):
	"""
		1. toMAF
		2. Combine Mafs per patient and assign scores.
		4. Combine All Mafs
	"""
	training_samples = {'TCGA-2H-A9GF', 'TCGA-2H-A9GO'}
	training_sample_list = [i for i in sample_list if i['PatientID'] in training_samples]
	prediction_sample_list = [i for i in sample_list if i['PatientID'] not in training_samples]
	#all_variants = GetVariantList(sample, options)
	#pprint(all_variants)

	truthsets = list()
	for training_case in training_samples:
		_subset_truthset = Truthset(training_case, options, training_type = training_type)
		truthsets.append(_subset_truthset)
	#Combine the truthsets


	SomaticSeqTraining(all_variants, options, truthset)
	SomaticSeqPredictor()


def SomaticSeqTraining(sample_variants, options, truthset):
	""" Trains the somaticseq algorithm.
		Parameters
		----------
			sample_variants: dict<>
				A dictionary mapping callers to the output vcfs.
		Inputs
		------

		Outputs
		-------------
			SOMATICSEQ_FOLDER/training
				Ensemble.sINDEL.tsv
				Ensemble.sINDEL.tsv.Classifier.RData
				Ensemble.sSNV.tsv
				Ensemble.sSNV.tsv.Classifier.RData

	"""
	variants = sample_variants
	output_folder = os.path.join(
		Pipeline_Folder, 'somaticseq', 'training-' + truthset.training_type)
	checkdir(output_folder)

	somaticseq_location = os.path.join(
		options['Programs']['SomaticSeq'],
		"SomaticSeq.Wrapper.sh")
	
	ada_script = os.path.join(
		options['Programs']['SomaticSeq'],
		'r_scripts',
		"ada_model_builder.R")
	
	gatk_location=options['Programs']['GATK']
	reference 	= options['Reference Files']['reference genome']
	cosmic 		= options['Reference Files']['cosmic']
	dbSNP 		= options['Reference Files']['dbSNP']

	print("Script Information: ")
	print("\t", ada_script)
	print("\t", os.path.isfile(ada_script))

	command = """{somaticseq} \
		--mutect2 {mutect} \
		--varscan-snv {varscan_snv} \
		--varscan-indel {varscan_indel} \
		--sniper {sniper} \
		--muse {muse} \
		--strelka-snv {strelka_snv} \
		--strelka-indel {strelka_indel} \
		--normal-bam {normal} \
		--tumor-bam {tumor} \
		--ada-r-script {ada} \
		--genome-reference {reference} \
		--cosmic {cosmic} \
		--dbsnp {dbSNP} \
		--gatk {gatk} \
		--inclusion-region {targets} \
		--truth-snv {truthset_snv} \
		--truth-indel {truthset_indel} \
		--output-dir {output_folder}""".format(
			somaticseq 		= somaticseq_location,
			gatk 			= gatk_location,
			ada 			= ada_script,

			normal 			= sample['NormalBAM'],
			tumor 			= sample['TumorBAM'],
			targets 		= sample['ExomeTargets'],
			reference 		= reference,
			cosmic 			= cosmic,
			dbSNP 			= dbSNP,

			muse 			= variants['muse'],
			mutect 			= variants['mutect'],
			sniper 			= variants['somaticsniper'],
			strelka_snv 	= variants['strelka-snv'],
			strelka_indel 	= variants['strelka-indel'],
			varscan_snv   	= variants['varscan-snv'],
			varscan_indel 	= variants['varscan-indel'],
			
			truthset_snv 	= truthset.snv_filename,
			truthset_indel 	= truthset.indel_filename,
			output_folder 	= output_folder)
	Terminal(command)

def SomaticSeqPredictor(sample_list, options):
	""" 
	"""
	somaticseq_output_folder = "/home/upmc/Documents/Genomic_Analysis/somaticseq"
	prediction_output_folder = os.path.join(somaticseq_output_folder, 'prediction')
	training_output_folder = os.path.join(somaticseq_output_folder, 'training')
	checkdir(prediction_output_folder)
	
	somaticseq_location = os.path.join(
		options['Programs']['SomaticSeq'],
		"SomaticSeq.Wrapper.sh")
	
	ada_script = os.path.join(options['Programs']['SomaticSeq'],
		"r_scripts",
		"ada_model_predictor.R")
	
	gatk_location = options['Programs']['GATK']

	normal = sample['NormalBAM']
	tumor = sample['TumorBAM']
	reference = options['Reference Files']['reference genome']
	cosmic = options['Reference Files']['cosmic']
	dbSNP = options['Reference Files']['dbSNP']

	snv_classifier = os.path.join(training_output_folder, "Ensemble.sSNV.tsv.Classifier.RData")
	indel_classifier = os.path.join(prediction_output_folder, "Ensemble.sINDEL.tsv.Classifier.RData")

	command = """{somaticseq} \
		--mutect2 {mutect} \
		--varscan-snv {varscan_snv} \
		--varscan-indel {varscan_indel} \
		--sniper {somaticsniper} \
		--muse {muse} \
		--strelka-snv {strelka_snv} \
		--strelka-indel {strelka_indel} \
		--normal-bam {normal} \
		--tumor-bam {tumor} \
		--ada-r-script {ada_model_predictor} \
		--classifier-snv {snv_classifier} \
		--classifier-indel {indel_classifier} \
		--pass-threshold 0.5 \
		--lowqual-threshold 0.1 \
		--genome-reference {reference} \
		--cosmic {cosmic} \
		--dbsnp {dbSNP} \
		--inclusion-region {targets} \
		--gatk {GATK} \
		--output-dir {output_folder}""".format(
			somaticseq = somaticseq_location,
			GATK = gatk_location,
			ada_model_predictor = ada_script,

			muse = variants['muse'],
			mutect = variants['mutect'],
			somaticsniper = variants['somaticsniper'],
			strelka_snv = variants['strelka-snv'],
			strelka_indel = variants['strelka-indel'],
			varscan_snv = variants['varscan-snv'],
			varscan_indel = variants['varscan-indel'],

			normal = normal,
			tumor = tumor,

			
			reference = reference,
			targets = sample['ExomeTargets'],
			cosmic = cosmic,
			dbSNP = dbSNP,
			snv_classifier = snv_classifier,
			indel_classifier = indel_classifier,
			output_folder = prediction_output_folder)
	Terminal(command)


options = configparser.ConfigParser()
options.read(OPTIONS_FILENAME)
"""
	Truthset -> Somaticseq training -> somaticseq prediction -> filter -> convert to MAF - > Combine MAFs
"""

class Pipeline:
	""" Processes the output of all callers into a final callset.
		Input
		-----
			sample files:
				SOMATIC_PIPELINE_FOLDER/barcode/caller/

		Output
		------
			Truthset: 
			SomaticSeq: 
		
		Overview
		--------
			Truthset -> Somaticseq training -> somaticseq prediction -> filter -> convert to MAF - > Combine MAFs

	"""
	def __init__(self, training_samples, prediction_samples, options, training_type):
		"""	Parameters
			----------
				training_samples: list<dict>
					A list of samples from the sample list to train SomaticSeq with.
				prediction_samples: list<dict>
					A list of samples to use with the SomaticSeq predictor.
				options: 
				training_type: {'Intersection', 'RNA'}
					The type of truthset to use.
		"""
		all_samples = training_samples + predition_samples
		#raw_variants = GetVariantList(sample, options, 'DNA-seq')

		truthset = Truthset(training_samples, options, training_type = training_type)

		SomaticseqTraining(training_samples, options, truthset)
		SomaticseqPredictor(prediction_samples, options)

		#Filter

		#Annotate and convert to MAF

		#Combine MAFs


if __name__ == "__main__" and True:
	with open(os.path.join(PIPELINE_FOLDER, "DNA-seq_Sample_List.tsv"), 'r') as file1:
		full_sample_list = list(csv.DictReader(file1, delimiter = '\t'))

	training_type = 'Intersection'
	training_samples = ['TCGA-2H-A9GF', 'TCGA-2H-A9GO']

	training_samples = [i for i in samples if i['PatientID'] in training_samples]
	prediction_samples = [i for i in samples if i['PatientID'] not in training_samples]

	Pipeline(training_samples, prediction_samples, options, training_type)

elif True:
	sample = {
		'PatientID':  "TCGA-2H-A9GF",
		'SampleID':   "TCGA-2H-A9GF-01A",
		'NormalID':   "TCGA-2H-A9GF-11A",
		'SampleUUID': "2d4f1ce4-4613-403a-90ec-fd6a551b6487",
		'NormalUUID': "2d1700f2-0f9b-4af6-a59b-4e6a51e01368",
		'NormalBAM':  "/home/upmc/Documents/genome_files/2d4f1ce4-4613-403a-90ec-fd6a551b6487/fcd15f0228c03bc602ae63c2e2b1b85c_gdc_realn.bam",
		'TumorBAM':   "/home/upmc/Documents/genome_files/2d1700f2-0f9b-4af6-a59b-4e6a51e01368/3e683408778f2ebe86f6c1f9ff936b0b_gdc_realn.bam",
		'ExomeTargets': "/home/upmc/Documents/Reference/SeqCap_EZ_Exome_v3_capture.hg38_GDC_official.bed"
	}
	#GATKCombineVariants(sample, options, variants, output_folder):
	output_folder = "/home/upmc/Documents/Genomic_Analysis/debug_folder/"
	CombineVariants(sample, options, output_folder)