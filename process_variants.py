import os
import collections
import sys

import configparser

import vcf
import csv
from pprint import pprint

"""
	1. Harmonize VCF fields
	2. Merge caller VCFs
	3. Convert VCFs to MAF
	4. Combine Patient MAFs
"""
GITHUB_FOLDER = os.path.join(os.getenv('HOME'), 'Documents', 'Github')
sys.path.append(GITHUB_FOLDER)
PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")

import pytools.systemtools as systemtools
import pytools.filetools as filetools
import pytools.tabletools as tabletools
import varianttools.callertools

class Workflow:
	def __init__(self, **kwargs):
		pass

class CombineCallsets:
	"""Combine separated indel and snv files
		WARNING: CombineVariants._copy_vcf() uses hard filter to remove variants with '/' formatting
	"""
	def __init__(self, sample, options, output_folder, variants):
		""" Parameters
			----------
				sample:
				options:
				output_folder:
				variants:
		"""
		self.cv_output_folder = output_folder
		self.gatk_program = options['Programs']['GATK']
		self.reference = options['Reference Files']['reference genome']

		modified_variants = self._modify_variants(variants)

		snv_variants, indel_variants = self.splitIndelSnvVariants(modified_variants) #saved in same folder as parent

		output_basename = os.path.join(output_folder, "{0}_vs_{1}.merged".format(sample['NormalID'], sample['SampleID']))
		self.snvs   = self.combineVariants(  snv_variants, output_basename + '.snp.vcf')
		self.indels = self.combineVariants(indel_variants, output_basename + '.indel.vcf')

	
	def _modify_merged_vcf(self, filename):
		basename = os.path.splitext(os.path.basename(filename))[0]
		basename = basename + ".modified.vcf"
		output_file = os.path.join(self.output_folder, basename)

		with open(filename, 'r') as vcf_file:
			reader = vcf.Reader(vcf_file)
			reader.infos['VAF'] = reader.formats['FREQ']._replace(type='Float')
			
			with open(output_file, 'w') as file2:
				writer = vcf.Writer(file2, reader)
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

		result = {
			'ALLELES': sample_alleles,
			'READS': sample_reads,
			'VAF': float("{0:.5f}".format(sample_vaf)),
		}
		return result

	def _modify_variants(self, vcfs):
		""" Some caller outputs are inconsistent and need to be modified.
		"""
		output_folder = self.cv_output_folder
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
				filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
				if filterOut:
					pass
				else:
					writer.write_record(record)
		return output_file
	def _modify_varscan_output(self, reader, filename):
		reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
		reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
		return reader
	def combineVariants(self, variants, output_file):
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

		order = "mutect,varscan,strelka,muse,somaticsniper" #ordered by VAF confidence
		
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
				gatk = 		self.gatk_program,
				reference = self.reference,
				muse = 		variants['muse'],
				mutect = 	variants['mutect'],
				ss = 		variants['somaticsniper'],
				strelka = 	variants['strelka'],
				varscan = 	variants['varscan'],
				rod = 		order,
				output = 	output_file)
		if not os.path.exists(output_file):
			systemtools.Terminal(command)
		return output_file
	def catVariants(self, left, right):
		""" Combines the SNV and Indel files. Assumes both are saved in the same folder. """
		l = os.path.splitext(os.path.splitext(left)[0])[0]
		output_file = l + '.cat.vcf'

		command = """java -cp {GATK} org.broadinstitute.gatk.tools.CatVariants \
		    -R {reference}\
		    -V {left} \
		    -V {right} \
		    -out {output}
		""".format(
			GATK = self.gatk_program,
			reference = self.reference,
			left = left,
			right = right,
			output = output_file)
		systemtools.Terminal(command)

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
					A list of samples to base the truthset on. Note: the sample variants will be found separately.
				options: dict<>
				training_type: {'Intersection', 'RNA-seq'}
				verbose: bool, default False
					Whether to print status messages.
		"""
		verbose = kwargs.get('verbose', False)
		kwargs['intersection'] = kwargs.get('intersection', 5)
		self.min_tumor_vaf = 0.08
		self.max_normal_vaf = 0.03

		self.truthset_folder = os.path.join(options['Pipeline Options']['processed vcf folder'], "truthset")
		self.gatk_program 	= options['Programs']['GATK']
		self.picard_program = options['Programs']['Picard']
		self.reference 		= options['Reference Files']['reference genome']


		
		self.training_type = training_type
		outputs = list()

		for sample in samples:
			if verbose:
				print(sample)
			if training_type == 'Intersection': sample_variants = GetVariantList(sample)
			elif training_type == 'RNA-seq': sample_variants = GetVariantList(sample, "RNA-seq")
			sample_truthset = self._per_sample(sample, options, training_type, **kwargs)
			outputs.append(sample_truthset)

		if len(outputs) == 1:
			sample_truthset = outputs[0]
		else:
			sample_truthset = self._combineTruthsets(outputs, options)
		#self.truthset = sample_truthset['truthset']
		self.filename       = sample_truthset['filename']
		self.indel_filename = sample_truthset['filename-indel']
		self.snv_filename   = sample_truthset['filename-snv']

	def __str__(self):
		string = "Truthset(type = {}, filename = {})".format(self.training_type, self.filename)
		return string
	def _per_sample(self, variants, options, training_type, **kwargs):
		""" Generates individual truthsets per sample.
			Available Keyword Arguments
			---------------------------
				'intersection': number of callers to count as the intersection.
				'variant_type': If not 'snv' or 'indel', the file will be split into
					an 'snv' and 'indel' file.
		"""
		normalID = sample['NormalID']
		tumorID = sample['SampleID']

		vcf_filenames = self._get_vcf_files(sample, options, training_type)

		output_vcf = os.path.join(                self.truthset_folder, "{0}_vs_{1}.{2}.truthset.vcf".format(normalID, tumorID, training_type))
		snv_filename = snv_vcf = os.path.join(self.truthset_folder, "{0}_vs_{1}.{2}.snv.truthset.vcf".format(normalID, tumorID, training_type))
		indel_filename = os.path.join(      self.truthset_folder, "{0}_vs_{1}.{2}.indel.truthset.vcf".format(NormalID, tumorID, training_type))
		output_folder = os.path.join(self.truthset_folder, 'merged_vcfs', sample['PatientID'])
		systemtools.checkDir(output_folder, True)
		
		combined_variants = CombineCallset(sample, options, output_folder, variants = variants)
		
		#combined_variants = CombineVariants(sample, options, output_folder)
		snv_variants = combined_variants.snvs
		indel_variants = combined_variants.indels
		
		#if not os.path.exists(output_vcf) or True:
		snv_truthset = self._generate_truthset(snv_variants, snv_filename, training_type, **kwargs)
		indel_truthset = self._generate_truthset(indel_variants, indel_filename, training_type, **kwargs)

		result = {
			'PatientID': sample['PatientID'],
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

		sampleIds = sorted(i['PatientID'] for i in samples)
		sampleIds = [i.split('-')[-1] for i in sampleIds] #Extracts patient id from the barcode. ex. TCGA-2H-A9GF -> A9GF
		sampleIds = ",".join(sampleIds)
		_folderName = lambda s: os.path.join(self.truthset_folder, s.format(sampleIds, self.training_type))

		filename       = _folderName("{0}.{1}.merged_truthset.vcf")
		indel_filename = _folderName("{0}.{1}.merged_truthset.indel.vcf")
		snv_filename   = _folderName("{0}.{1}.merged_truthset.snv.vcf")

		snv_variants   = "--variant " + " --variant ".join(i['filename-snv'] for i in samples) + ' \\'
		indel_variants = "--variant " + " --variant ".join([i['filename-indel'] for i in samples]) + '\\'
		snv_variants   = " ".join(["I={0}".format(i['filename-snv']) for i in samples])
		indel_variants   = " ".join(["I={0}".format(i['filename-indel']) for i in samples])
		
		gatk_base_command = """ java -cp {program} org.broadinstitute.gatk.tools.CatVariants \
				    -R {reference} \
				    {variants}
				    -out {output}"""
		
		picard_base_command = """java -jar {program} SortVcf {variants} O={output}"""
		snv_command = picard_base_command.format(
			program 	= self.picard_program,
			reference 	= self.reference,
			variants 	= snv_variants,
			output 		= snv_filename)

		indel_command = picard_base_command.format(
			program 	= self.picard_program,
			reference 	= self.reference,
			variants 	= indel_variants,
			output 		= indel_filename)

		if not os.path.exists(snv_filename):
			systemtools.Terminal(snv_command)
		if not os.path.exists(indel_filename):
			systemtools.Terminal(indel_command)
		result = {
			'PatientID': sampleIds,
			'filename': filename,
			'filename-indel': indel_filename,
			'filename-snv': snv_filename
		}

		return result
	
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
		_somatic_vaf = (tumor_vaf >= self.min_tumor_vaf and normal_vaf < self.max_normal_vaf) or (tumor_vaf < self.min_tumor_vaf and normal_vaf == 0.0)
		if _somatic_vaf and not filter_out:
			validation_status = 1 #Somatic
		elif (tumor_vaf < self.min_tumor_vaf and (normal_vaf > 0.0 and normal_vaf < self.max_normal_vaf)) or filter_out:
			validation_status = 0 #Non-Somatic
		else: validation_status = 2 #'Unknown'

		validation_status = int(validation_status == 1)

		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'VAF',
			'validation status': validation_status
		}
		return row
	

class VCFtoMAF:
	def __init__(self, sample, input_vcf, options, truthset = None):
		""" Converts a VCF file to and annotated MAF file. Include additional fields such as VAF.
			Parameters
			----------
				sample:
				input_vcf: string, dict<string:string>
					The vcf file to convert.
				options:
				truthset:
			Additional Columns
			------------------
				VAF: The variant allele frequency

		"""
		self.reference_build 	= 'GRCh38'
		self.vcftomaf_program 	= options['Programs']['vcf2maf']
		self.vep_program 		= options['Programs']['varianteffectpredictor']
		self.reference 			= options['Reference Files']['reference genome']
		self.maf_fieldnames 	= [ 
			"Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", 
			"End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
			"Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", 
			"Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
			"Tumor_Validation_Allele1", "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
			"Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status",
			"Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", "BAM_File", "Sequencer",
			"Tumor_Sample_UUID", "Matched_Norm_Sample_UUID", "HGVSc", "HGVSp", "HGVSp_Short", "Transcript_ID",
			"Exon_Number", "t_depth", "t_ref_count", "t_alt_count",  "n_depth", "n_ref_count", "n_alt_count",
			"all_effects", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", 
			"CDS_position",	"Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM"
			"DISTANCE", "STRAND_VEP", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "CANONICAL", "CCDS",
			"ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "RefSeq", "SIFT", "PolyPhen", "EXON", "INTRON",
			"DOMAINS", "GMAF", "AFR_MAF", "AMR_MAF", "ASN_MAF", "EAS_MAF", "EUR_MAF", "SAS_MAF", "AA_MAF",
			"EA_MAF", "CLIN_SIG", "SOMATIC", "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS",
			"MOTIF_SCORE_CHANGE", "IMPACT", "PICK", "VARIANT_CLASS", "TSL", "HGVS_OFFSET", "PHENO",
			"MINIMISED", "ExAC_AF", "ExAC_AF_AFR", "ExAC_AF_AMR", "ExAC_AF_EAS", "ExAC_AF_FIN", 
			"ExAC_AF_NFE", "ExAC_AF_OTH", "ExAC_AF_SAS", "GENE_PHENO", "FILTER", "flanking_bps", "variant_id",
			"variant_qual", "ExAC_AF_Adj", "ExAC_AC_AN_Adj", "ExAC_AC_AN", "ExAC_AC_AN_AFR", "ExAC_AC_AN_AMR",
			"ExAC_AC_AN_EAS", "ExAC_AC_AN_FIN", "ExAC_AC_AN_NFE", "ExAC_AC_AN_OTH", "ExAC_AC_AN_SAS", "ExAC_FILTER"
		]

		self.output_folder = options['Pipeline Options']['maf folder']
		raw_maf = os.path.join(output_folder, "{0}.raw.maf".format(sample['PatientID']))
		self.maf = os.path.join(output_folder, "{0}.maf".format(sample['PatientID']))
		maf_file = self.vcftomaf(sample, input_vcf, raw_maf)
		#maf_file = self.modifyMAF(sample, options, maf_file, truthset)

		print("Output MAF: ", self.final_maf)


		#file1.write("#version 2.4\n")
		#writeTSV(maf_file, output_filename, self.maf_fieldnames)
		#self.filename = output_filename

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
	def vcftomaf(sample, input_vcf, output_maf):
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
				script 		= self.vcftomaf_program,
				vep 		= self.vep_program,
				reference 	= self.reference,
				ncbi 		= self.reference_build,
				tumor 		= sample['SampleID'],
				normal 		= sample['NormalID'],
				vcf 		= input_vcf,
				maf 		= output_maf)
		#print(command)
		if not os.path.isfile(output_maf):
			systemtools.Terminal(command)

		return output_maf

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


class SomaticSeq:
	def __init__(self, training_samples, prediction_samples, options, truthset):
		"""
			Generate tables for the training samples individually.
			Combine the tables for the training samples.
			Train the model using the combined output generated previously.
			Predict varaints for all prediction samples.

			Parameters
			----------
				training_samples: list<dict>
				prediction_samples: list<dict>
				options:
				truthset: string, Truthset
					Either a training type to generate a truthset with or a pre-configured truthset.

			Trainer Output
			--------------
				[training folder]
					files
					[tables]
						[patient barcode]
							files
				[prediction folder]
					[patient barcode]
						files
		"""

		self.somaticseq_folder = options['Programs']['SomaticSeq']
		self.somaticseq_program = os.path.join(self.somaticseq_folder, "SomaticSeq.Wrapper.R")
		self.ada_builder_script = os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_builder.RData")
		self.ada_predictor_script=os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_predictor.RData")
		self.tsv_to_vcf_script  = os.path.join(self.somaticseq_folder, "SSeq_tsv2vcf.py")
		self.gatk_program = options['Programs']['GATK']

		self.cosmic 	= options['Reference Files']['cosmic']
		self.dbSNP 		= options['Reference Files']['dbSNP']
		self.reference 	= options['Reference Files']['reference genome']

		self.somaticseq_output_folder = os.path.join(PIPELINE_FOLDER, 'somaticseq')

		#Generate the truthset.
		if isinstance(truthset, str):
			#The truthset training type was supplied instead of the truthset itself.
			self.training_type = truthset
			self.truthset = Truthset(training_samples, options, training_type)
		else:
			#A truthset object was passed
			self.training_type = truthset.training_type
			self.truthset = truthset
		print(Truthset)

		self.training_folder   = os.path.join(self.somaticseq_output_folder, 'training-'   + self.training_type)
		self.prediction_folder = os.path.join(self.somaticseq_output_folder, 'prediction-' + self.training_type)
		filetools.checkDir(self.training_folder)
		filetools.checkDir(self.prediction_folder)


		#Generate separate tables for each training sample.
		if   self.training_type == 'Intersection': file_type = 'DNA-seq'
		elif self.training_type == 'RNA-seq':      file_type = 'RNA-seq'

		self.trainers = list()
		for training_sample in training_samples:
			sample_variants = GetVariantList(training_sample, file_type)
			trainer = self._runTrainer(training_sample, sample_variants)
			self.trainers.append(trainer)

		#self.ensemble_tables = self._combineTrainers(self.trainers)
		self.ensemble_tables = self._combineTables(self.trainers)
		self.trainer = self._runManualClassifier(
			self.ensemble_tables['indelTable'], 
			self.ensemble_tables['snvTable'])

		#Run the prediction model on the prediction variants.
		self.predictions = list()
		for prediction_sample in prediction_samples:
			sample_variants = GetVariantList(prediction_sample)
			prediction = self._runPredictor(prediction_sample, trainer)
			self.predictions.append(prediction)


	def _runTrainer(self, sample, variants):
		"""
			Parameters
			----------
				sample: dict<>
					The sample's row from the sample list. Contains information related to the genome files.
				variants: dict<>
					Each callset generated for the sample.
				truthset: Truthset
					The truthset for the analysis. Should have both SNV and Indel values.
			Outputs
			-------
				[somaticseq output folder]/training-[training type]/[barcode]/[files]
		"""
		output_folder = os.path.join(self.training_folder, 'tables', sample['PatientID'])
		filetools.checkDir(output_folder, True)
		#			--ada-r-script {ada} \
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
			--genome-reference {reference} \
			--cosmic {cosmic} \
			--dbsnp {dbSNP} \
			--gatk {gatk} \
			--inclusion-region {targets} \
			--truth-snv {truthset_snv} \
			--truth-indel {truthset_indel} \
			--output-dir {output_folder}""".format(
				somaticseq 		= self.somaticseq_wrapper,
				gatk 			= self.gatk_program,
				ada 			= self.ada_builder_script,

				normal 			= sample['NormalBAM'],
				tumor 			= sample['TumorBAM'],
				targets 		= sample['ExomeTargets'],
				reference 		= self.reference,
				cosmic 			= self.cosmic,
				dbSNP 			= self.dbSNP,

				muse 			= variants['muse'],
				mutect 			= variants['mutect'],
				sniper 			= variants['somaticsniper'],
				strelka_snv 	= variants['strelka-snv'],
				strelka_indel 	= variants['strelka-indel'],
				varscan_snv   	= variants['varscan-snv'],
				varscan_indel 	= variants['varscan-indel'],
				
				truthset_snv 	= self.truthset.snv_filename,
				truthset_indel 	= self.truthset.indel_filename,
				output_folder 	= output_folder)

		expected_output = {
			'indelTable': os.path.join(output_folder, "Ensemble.sINDEL.tsv"),
			'snvTable': os.path.join(ouput_folder, "Ensemble.sSNV.tsv"),
			'indelClassifier': os.path.join(output_folder, "Ensemble.sINDEL.tsv.Classifier.RData"),
			'snvClassifier': os.path.join(ouput_folder, "Ensemble.sSNV.tsv.Classifier.RData")
		}
		if any([not os.path.exists(fn) for fn in expected_output.values()]):
			systemtools.Terminal(command)
		expected_output['PatientID'] = sample['PatientID']
		return expected_output

	def _TsvToVcf(self, filename):
		""" Converts a SomaticSeq TSV file to a VCF file in the same folder.
		"""
		output_file = os.path.splitext(vcf_filename)[0] + '.vcf'

		command = """{script} -tsv {table} -vcf {output} -pass 0.7 -low 0.1 -tools {tools} -phred"""
		command = command.format(
			script = self.tsv_to_vcf_script,
			table = filename,
			output = output_file,
			tools = "")

		if not os.path.exists(output_file):
			systemtools.Terminal(command)

		return output_file

	def _combineTables(self, tables):
		""" Combines the output Ensemble tables form somaticseq.
		"""
		output_folder = os.path.join(self.training_folder, "combined_tables")
		filetools.checkDir(output_folder)
		training_ids = sorted(i['PatientID'] for i in tables)
		training_ids_string = ','.join([i.split('-')[2] for i in training_ids])
		indel_table_filename = os.path.join(output_folder, "Ensemble.sINDEL.combined.{0}.tsv".format(training_ids_string))
		snv_table_filename = os.path.join(output_folder, "Ensemble.sSNV.combined.tsv.{0}.tsv".format(training_ids_string))

		indel_table = list()
		snv_table = list()
		pprint(tables)
		for filename in [i['indelTable'] for i in tables]:
			table, indel_fieldnames = tabletools.readCSV(filename, True)
			indel_table += table
		for filename in [i['snvTable'] for i in tables]:
			table, snv_fieldnames = tabletools.readCSV(filename, True)
			snv_tables += table

		tabletools.writeCSV(indel_table, indel_table_filename, fieldnames = indel_fieldnames)
		tabletools.writeCSV(snv_table, snv_table_filename, fieldnames = snv_fieldnames)

		output = {
			'snvTable': snv_table_filename,
			'indelTable': indel_table_filename
		}
		return output

	def _runManualClassifier(self, indel_table, snv_table):
		""" Trains the model using the provided tables.
		"""
		output_folder = os.path.dirname(snv_table)
		command = "{script} {table}"
		indel_command = command.format(self.ada_builder_script, indel_table)
		snv_command = command.format(self.ada_builder_script, snv_table)

		expected_output = {
			'indelTable': os.path.join(output_folder),
			'snvTable': os.path.join(output_folder),
			'indelClassifer': os.path.join(output_folder),
			'snvClassifier': os.path.join(output_folder)
		}

		if not os.path.exists(expected_output['indelClassifier']):
			os.system(indel_command)
		if not os.path.exists(expected_output['snvClassifer']):
			os.system(snv_command)

		return expected_output

	def _runPredictor(self, sample, trainer):
		""" Runs the Somaticseq prediction model.
		"""
		output_folder = os.path.join(self.prediction_folder, sample['PatientID'])
		normal = sample['NormalBAM']
		tumor = sample['TumorBAM']

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
				somaticseq 			= self.somaticseq_program,
				GATK 				= self.gatk_program,
				ada_model_predictor = self.ada_predictor_script,

				muse 		  = variants['muse'],
				mutect 		  = variants['mutect'],
				somaticsniper = variants['somaticsniper'],
				strelka_snv   = variants['strelka-snv'],
				strelka_indel = variants['strelka-indel'],
				varscan_snv   = variants['varscan-snv'],
				varscan_indel = variants['varscan-indel'],

				normal = normal,
				tumor = tumor,

				
				reference 	= self.reference,
				targets 	= sample['ExomeTargets'],
				cosmic 		= self.cosmic,
				dbSNP 		= self.dbSNP,

				snv_classifier   = trainer['snvClassifier'],
				indel_classifier = trainer['indelClassifier'],
				output_folder = self.prediction_folder)
		indel_table = os.path.join(output_folder, "Ensemble.sINDEL.tsv")
		snv_table = os.path.join(output_folder, "Ensemble.sSNV.tsv")
		if not os.path.exists(indelTable) or not os.path.exists(snvTable):
			systemtools.Terminal(command)

		indel_vcf =	self.TsvToVcf(indelTable)
		snv_vcf = self.TsvToVcf(snvTable)

		expected_output = {
			'PatientID': sample['PatientID'],
			'indelTable': indel_table,
			'snvTable': snv_table,
			'indelVCF': indel_vcf,
			'snvVCF': snv_vcf,
			'indelClassifer': os.path.join(output_folder),
			'snvClassifier': os.path.join(output_folder)
		}

		return expected_output

class FilterVCF:
	"""
		Filters variants based on a set of common filters.
		Generates two files, one containing the variants that passed all filters and one with variants
		that were rejected.

		Filters
		-------
			a. filter variants present in dbSNP, 1000 genomes project, etc.
			b. filter out known indels.
			c. filter out germline variants.
			d. filter out variants where the base quality is low
			e. filter out variants based on mutation consequence.
			f. Annotate variants
			g. somatic coverage filtering.
	"""
	def __init__(self, vcf_file, options, isIndel):
		"""
			vcf_file: Path to the vcf file to filter.
			options:
			isIndel: Needed to use indel- or snv- specific filters.
		"""
		self.dbSNP = options['Reference Files']['dbSNP']
		self.output = vcf_file

	def _filterVcfFile(self, filename):
		output_folder = ""
		passed_output = ""
		rejected_output=""

		num_passed_records = 0
		num_rejected_records=0

		passed_records = list()
		rejected_records= list()
		with open(filename, 'r') as vcf_file:
			reader = vcf.Reader(vcf_file)
			for record in reader:
				filters = self._getFilterStatus(record)
				if filters['filterOut']:
					num_rejected_records += 1
					rejected_records.append(record)
				else:
					num_passed_records += 1
					passed_records.append(record)
			self.writeVCF(passed_records, passed_output, reader)
			self.writeVCF(rejected_records, rejected_output, reader)

		result = {
			'numPassedRecords': num_passed_records,
			'numRejectedRecords': num_rejected_records,
			'passedVCF': passedOutput,
			'rejectedVCF': rejectedOutput
		}

		return result
	def writeVCF(self, records, filename, reader):
		with open(filename, 'w') as file1:
			writer = vcf.Writer(file1, reader)
			for record in records:
				writer.write_record(record)
		return filename

	def _getFilterStatus(self, record, record_type):
		"""
			a. filter variants present in dbSNP, 1000 genomes project, etc.
			b. filter out known indels.
			c. filter out germline variants.
			d. filter out variants where the base quality is low
			e. filter out variants based on mutation consequence.
			f. Annotate variants
			g. somatic coverage filtering.
		"""
		dbSNP_status       = self._filterFromdbSNP(record)
		known_indel_status = self._filterKnownIndels(record)
		base_quality_status= self._filterBaseQuality(record)



	def _filterFromdbSNP(self, record):
		pass

	def _filterKnownIndels(self, record):
		pass

	def _filterBaseQuality(self, record):
		pass

	def _filterConsequence(self, record):
		pass

	def _filterGermlineVariants(self, record):
		pass


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
			Folder Setup:
				1_input_files
					Case1
						Callset1.vcf
						Callset2.vcf
						...
					Case2
					Case3
					...
				2_processed_vcfs
					truthset
						truthset1
						truthset2
						...
				somaticseq
					training
						..files
					prediction
						..files



	"""
	def __init__(self, training_samples, prediction_samples, options, training_type, **kwargs):
		"""	Parameters
			----------
				training_samples: list<dict>
					A list of samples from the sample list to train SomaticSeq with.
				prediction_samples: list<dict>
					A list of samples to use with the SomaticSeq predictor.
				options: config
				training_type: {'Intersection', 'RNA'}
					The type of truthset to use.
						'Intersection': Will use the intersection of the callsets
						'RNA': Will use variants called via RNA-seq as the truthset.
			Keyword Arguments
			-----------------
				verbose: bool; default True
					If true, prints status messages.
		"""
		verbose = kwargs.get('verbose', True)

		if verbose:
			print("training samples: ", ', '.join(i['PatientID'] for i in training_samples))
			print("prediction samples: ", ', '.join(i['PatientID'] for i in prediction_samples))

		all_samples = training_samples + prediction_samples
		
		#raw_variants = GetVariantList(sample, options, 'DNA-seq')

		truthset = Truthset(training_samples, options, training_type = training_type, **kwargs)

		somaticseq = SomaticSeq(training_samples, prediction_samples, options, truthset = truthset)
		predictions = somaticseq.predictions #list of dictionaries with the output filenames for a given sample.
		
		#Filter
		filtered_samples = list()
		for sample_vcfs in predictions:
			sample = [i for i in prediction_samples if i['PatientID'] == sample_vcfs['PatientID']][0]
			indel_vcf = sample_vcfs['indelVCF']
			snv_vcf   = sample_vcfs['snvVCF']
			sample_vcfs['indelFilteredVCF'] = FilterVCFs(indel_vcf, options).output
			sample_vcfs['snvFilteredVCF']   = FilterVCFs(snv_vcf, options).output

			filtered_samples.append(sample_vcfs)
		
		#Annotate and convert to MAF
		maf_samples = list()
		for filtered_sample in filtered_samples:
			sample = [i for i in prediction_samples if i['PatientID'] == filtered_sample['PatientID']]
			indel_vcf = filtered_sample['indelFilteredVCF']
			snv_vcf = filtered_sample['snvFilteredVCF']
			indel_maf = VCFtoMAF(sample, indel_vcf, optionss)
			snv_maf = VCFtoMAF(sample, snv_vcf, options)
			filtered_sample['indelMAF'] = indel_maf
			filtered_sample['snvMAF'] = snv_maf
			maf_samples.append(filtered_sample)
		
		#Combine MAFs




"""
	Truthset -> Somaticseq training -> somaticseq prediction -> filter -> convert to MAF - > Combine MAFs
"""

if __name__ == "__main__" and True:
	options = configparser.ConfigParser()
	options.read(OPTIONS_FILENAME)

	documents_folder = os.path.join(os.getenv('HOME'), 'Documents')

	full_sample_list_filename = os.path.join(documents_folder, "DNA-seq_Sample_List.tsv")
	full_sample_list = tabletools.readCSV(full_sample_list_filename)

	training_type = 'Intersection'
	training_samples = ['TCGA-2H-A9GF', 'TCGA-2H-A9GO']


	training_samples 	= [i for i in full_sample_list if i['PatientID'] in training_samples]
	prediction_samples 	= [i for i in full_sample_list if i['PatientID'] not in training_samples]

	Pipeline(training_samples, prediction_samples, options, training_type)

elif False:
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
	CombineCallset(sample, options, output_folder)
	#1940