import os
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
if os.name == 'nt':
	GITHUB_FOLDER = os.path.join(os.getenv('USERPROFILE'), 'Documents', 'Github')
else:
	GITHUB_FOLDER = os.path.join(os.getenv('HOME'), 'Documents', 'Github')
sys.path.append(GITHUB_FOLDER)

from pipeline_manager import getPipelineFolder, getSampleCallset
from somaticseq import SomaticSeqPipeline

import pytools.systemtools as systemtools

import pytools.tabletools as tabletools
import varianttools.callertools as callertools
import varianttools.vcftools as vcftools


PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")

# For testing
# set to None to disable
TEST_CHROMOSOME = 'chr1'
# Relevant Folders
VCF_FOLDER = os.path.join(PIPELINE_FOLDER, "1_vcf_files")
TRUTHSET_FOLDER = os.path.join(PIPELINE_FOLDER, "2_truthsets")



class GATKMergeSampleCallsets:
	"""Merges callsets via GATK CombineVariants.
		WARNING: CombineVariants._copy_vcf() uses hard filter to remove variants with '/' formatting
		Usage
		-----
			merger = MergeSampleCallsets()
			output = merger(sample, output_folder, callset)
		Requirements
		------------
			GATK
			Reference Genome
	"""
	def __init__(self, merge_options, **kwargs):
		""" Parameters
			----------
				sample:
				merge_options:
				output_folder:
				variants:
		"""
		log_message = "GATKMergeSampleCallsets(merge_options = {})".format(merge_options)
		print(log_message)
		print('\tKeyword Arguements:')
		for k, v in kwargs.items():
			print("\t{}\t{}".format(k, v))
		if merge_options is not None:
			self.gatk_program = merge_options['Programs']['GATK']
			self.reference = merge_options['Reference Files']['reference genome']
		else:
			self.gatk_program = kwargs.get('GATK')
			self.reference = kwargs.get('reference')
		if self.gatk_program is None or self.reference is None:
			message = "GATK only runs on linux!"
			raise OSError(message)

	def __call__(self, callset, **kwargs):
		""" Merges the provided callset
			Keyword Arguments
			-----------------
				* 'filename': The filename to save the merged callset as.
				* 'kind': If filename is not provided, 'kind' indicates the patient
					folder to use. Both 'patientId' and 'kind' must be provided.

		"""
		output_filename = kwargs.get('filename')
		self._checkCallsetFormat(callset)
		output_filename = self.gatkCombineVariants(callset, output_filename)
		return output_filename
			
	def _modify_merged_vcf(self, filename):
		""" Adds the VAF to the 'Info' field of the output file.
			The new file will be saved to the same folder as the original.
		:param filename: string
			Path to the merged file.
		:return: string
			path to the output file.
		"""
		output_folder = os.path.dirname(filename)
		basename = os.path.splitext(os.path.basename(filename))[0]
		basename = basename + ".modified.vcf"
		output_file = os.path.join(output_folder, basename)

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

	@staticmethod
	def _checkCallsetFormat(callset):
		""" Ensures a callset is formatted as a dictionary with a single file
			per caller. 
		"""
		caller_name_format_status = any('-' in c for c in callset.keys())

		callset_valid = not caller_name_format_status
		if not callset_valid:
			message = "The callset names include '-'"
			raise ValueError(message)
		return caller_name_format_status

	@staticmethod
	def _getVAF(record):
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
		record_sample = [i for i in record.samples if i.sample == 'TUMOR'].pop()
		if 'FREQ' in info_keys:
			VAF = float(record_sample['FREQ'][:-1])
		elif 'DP' in info_keys and 'AD' in info_keys:
			reads = record_sample['DP']
			alleles = record_sample['AD'][1]
			VAF = (alleles / reads) * 100
		elif 'AF' in info_keys:
			VAF = record_sample['AF']
		elif 'DP4' in info_keys:
			reads = sum(record_sample['DP4'])
			alleles = sum(record_sample['DP4'][2:])
			VAF = alleles / reads
		else:
			alleles = [i for i in ['A', 'C', 'G', 'T'] if i != record.REF]
			sample_reads = sum([record_sample[i+'U'][1] for i in (alleles + [record.REF])])
			sample_alleles = sum([record_sample[i+'U'][1] for i in alleles])
			if sample_reads == 0: sample_vaf = 0
			else:
				sample_vaf = sample_alleles / sample_reads
			VAF = sample_vaf
		return VAF

	def _modify_variants(self, vcfs):
		""" Some caller outputs are inconsistent and need to be modified.
		"""
		output_variants = dict()
		for caller, vcf_filename in vcfs.items():
			current_output_folder = os.path.dirname(vcf_filename)
			basename = os.path.splitext(os.path.basename(vcf_filename))[0] + ".modified.vcf"
			output_file = os.path.join(current_output_folder, basename)

			with open(vcf_filename, 'r') as vcf_file:
				reader = vcf.Reader(vcf_file)
				if 'varscan' in caller:
					reader = self._modify_varscan_output(reader)
				output_variants[caller] = self._copy_vcf(reader, output_file)

		return output_variants
	@staticmethod
	def _copy_vcf(reader, output_file):
		with open(output_file, 'w') as file1:
			writer = vcf.Writer(file1, reader)
			for record in reader:
				filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
				if filterOut:
					pass
				else:
					writer.write_record(record)
		return output_file

	@staticmethod
	def _modify_varscan_output(reader):
		reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
		reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
		return reader

	def _combineSplitVariants(self, patientId, output_folder, callset):
		""" Merges the indel/snp variants separately.
			Parameters
			----------
				sample: dict<>, string
					sample information as a dict or the patient id.
				output_folder: path
					Folder to save the merged variant files in.
				callset: Callset; default None
					If not provided, a generic callset will be created.
		"""
		if isinstance(patientId, dict):
			patientId = patientId['PatientID']

		snp_variants   = callset('snp')
		indel_variants = callset('indel')

		basename = "{}.merged_variants".format(patientId)
		output_snp_filename   = os.path.join(output_folder, basename + '.snp.vcf')
		output_indel_filename = os.path.join(output_folder, basename + '.indel.vcf')

		merged_snp_filname    = self.gatkCombineVariants(snp_variants,   output_snp_filename)
		merged_indel_filename = self.gatkCombineVariants(indel_variants, output_indel_filename)

		result = {
			'snp': merged_snp_filname,
			'indel': merged_indel_filename
		}
		return result

	def gatkCombineVariants(self, variants, output_file):
		""" Uses GATK CombineVariants to merge the calls from each caller into a single file.
			Parameters
			----------
				variants: dict<caller, path>
					A dictionary linkng each caller to its harmonized output.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Returns
			-------
				Output_file: string
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.merged.vcf
		"""

		order = "mutect2,varscan,strelka,muse,somaticsniper"  # ordered by VAF confidence
		variant_command = ['--variant:{} "{}"'.format(k, v) for k, v in variants.items()]
		variant_command = ' \\\n'.join(variant_command)
		command = """java -jar "{gatk}" \
			-T CombineVariants \
			-R "{reference}" \
			{variants} \
			-o "{output}" \
			-genotypeMergeOptions PRIORITIZE \
			-priority {rod}"""
		command = command.format(
				gatk =      self.gatk_program,
				reference = self.reference,
				variants =  variant_command,
				rod =       order,
				output =    output_file)

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
	def __init__(self, samples, truthset_options, truthset_type, **kwargs):
		""" Generates a truthset based on one or more samples.
			Option 1: intersection of all 5 callers.
			Option 2: Dream-SEQ data
			Option 3: RNA-seq?
			Parameters
			----------
				samples: dict<>
					A list of samples to base the truthset on. 
					The sample variants will be found separately.
					Only the original callset can be used to generate the truthset.
				truthset_options: dict<>
				truthset_type: {'Intersection', 'RNA-seq'}
				verbose: bool, default False
					Whether to print status messages.
			Keyword Arguments
			-----------------
				'indel_intersection': int; default 2
					The number of callers that must have called an indel
					variant to have it marked as a true positive.
				'snp_intersection': int; default 5
					The number of callers that must have called a snp
					variant to have it marked as a true positive.
		"""
		self.debug = True
		if self.debug:

			print("\ttruthset_type = {}".format(truthset_type))

		# Parameters to use when generating the truthsets.
		self.truthset_type = truthset_type

		self.indel_intersection = kwargs.get('indel_intersection', 2)
		self.snp_intersection   = kwargs.get(  'snp_intersection', 5)

		self.min_tumor_vaf  = 0.08
		self.max_normal_vaf = 0.03

		self.gatk_program   = truthset_options['Programs']['GATK']
		self.picard_program = truthset_options['Programs']['Picard']
		self.reference      = truthset_options['Reference Files']['reference genome']
		
		outputs = list()
	
		# Generate a separate truthset for each sample, then merge.
		if not isinstance(samples, list): samples = [samples]
		for sample in samples:
			sample_truthset = self._generateSampleTruthset(sample, **kwargs)
			outputs.append(sample_truthset)

		if len(outputs) == 1:
			sample_truthset = outputs[0]
		else:
			sample_truthset = self._combineTruthsets(outputs)

		self.indel_filename = sample_truthset['filename-indel']
		self.snp_filename   = sample_truthset['filename-snp']

	def _generateSampleTruthset(self, sample, **kwargs):
		""" Generates individual truthsets per sample.
			Available Keyword Arguments
			---------------------------
				'intersection': number of callers to count as the intersection.
				'variant_type': If not 'snv' or 'indel', the file will be split into
					an 'snv' and 'indel' file.
		"""
		patientId = sample['PatientID']
		################################# Define Filenames ########################################
		base_truthset_files_folder = getPipelineFolder(step = 'truthset', patientId = patientId)

		merged_indel_callset_filename = os.path.join(
			base_truthset_files_folder,
			"{}.raw.merged.indel.vcf".format(patientId)
		)
		merged_snp_callset_filename = os.path.join(
			base_truthset_files_folder,
			"{}.raw.merged.snp.vcf".format(patientId)
		)
		
		if self.truthset_type == 'Intersection':
			snp_truthset_label = "intersectionOf{}Callers".format(self.snp_intersection)
			indel_truthset_label="intersectionOf{}Callers".format(self.indel_intersection)
		else:
			message = "Truthset type wasn't implemented."
			raise NotImplementedError(message)


		final_indel_truthset_filename = os.path.join(
			base_truthset_files_folder,
			"{}.{}.final.indel.{}.truthset.vcf".format(patientId, training_type, indel_truthset_label)
		)

		final_snp_truthset_filename = os.path.join(
			base_truthset_files_folder,
			"{}.{}.final.snp.{}.truthset.vcf".format(patientId, training_type, snp_truthset_label)
		)

		###### Merge the callset files. The final truthset will be generated from this file. ######
		# Retrieve the callset for this sample.
		# Want to use the files that have been split into indels/snps and corrected.

		sample_indel_callset = getSampleCallset(patientId, 'original-fixed-split', 'indel')
		sample_snp_callset   = getSampleCallset(patientId, 'original-fixed-split', 'snp')
		
		if self.debug:
			print("\tTruthset._getSampleCallset: number of snp files found: ", len(sample_snp_callset))
			print("\tTruthset._getSampleCallset: number of indel files found: ", len(sample_indel_callset))
		
		if len(sample_snp_callset) < 5 or len(sample_indel_callset) < 5:
			pprint(sample_indel_callset)
			pprint(sample_snp_callset)
			raise ValueError()

		# Merge the full outputs from each caller.
		# The snps and indels need to be merged separately of each other, due to varscan/strelka
		merged_indel_callset_filename = GATK_MERGE_CALLSET(
			callset = sample_indel_callset,
			filename = merged_indel_callset_filename
		)
		merged_snp_callset_filename   = GATK_MERGE_CALLSET(
			callset = sample_snp_callset,
			filename = merged_snp_callset_filename
		)
		
		# Generate two truthsets, one for snps and one for indels.
		final_indel_truthset_filename = self._generate_truthset(
			input_vcf = merged_indel_callset_filename,
			output_vcf = final_indel_truthset_filename,
			training_type = self.truthset_type,
			**kwargs
		)

		final_snp_truthset_filename = self._generate_truthset(
			input_vcf = merged_snp_callset_filename,
			output_vcf = final_snp_truthset_filename,
			training_type = self.truthset_type,
			**kwargs
		)

		result = {
			'patientId':      sample['PatientID'],
			'filename-indel': final_indel_truthset_filename,
			'filename-snp':   final_snp_truthset_filename
		}
		return result

	@staticmethod
	def _getSampleVAF(record):
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

		for sample in record.samples:
			sample_fields = sample.data._asdict()

			if 'DP4' in sample_fields:          caller = 'somaticsniper'
			elif 'ALT_F1R2' in sample_fields:   caller = 'mutect'
			elif 'FREQ' in sample_fields:       caller = 'varscan'
			elif 'DP' in sample_fields:         caller = 'muse'
			else:                               caller = 'strelka'
			caller = caller.lower()

			if caller == 'somaticsniper':
				
				sample_reads = sum(sample_fields['DP4'])
				sample_alleles = sum(sample_fields['DP4'][2:])
				sample_vaf = sample_alleles / sample_reads

			elif caller == 'mutect':
				sample_reads =  sample_fields['ALT_F1R2'] + sample_fields['ALT_F2R1']
				sample_reads += sample_fields['REF_F1R2'] + sample_fields['REF_F2R1']
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
				message = "WARNING (getSampleVAF): {0} is not a caller!".format(caller)
				raise ValueError(message)

			result = {
				'reads': sample_reads,
				'alleles': sample_alleles,
				'vaf': sample_vaf
			}
			response[sample.sample] = result

		return response

	def _generate_truthset(self, input_vcf, output_vcf, **kwargs):
		""" Generates a truthset
			Parameters
			----------
				sample
				options
				training_type: {'RNA-seq', 'VAF', 'Intersection'}
		 """
		if self.debug:
			print("Generating Truthset...")
			print("\tInput: ", input_vcf)
			print("\tOutput:", output_vcf)
		with open(input_vcf, 'r') as input_file:
			reader = vcf.Reader(input_file)
			with open(output_vcf, 'w') as output_file:
				writer = vcf.Writer(output_file, reader)

				for index, record in enumerate(reader):
					recordStatus = self._getRecordValidationStatus(record, **kwargs)
					isValid = recordStatus['validation status']
					
					if isValid:
						writer.write_record(record)
		return output_vcf

	def _getRecordValidationStatus(self, record, **kwargs):
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
			# Assume the record is from the vcf from the RNA-seq pipeline.
			# Assume all RNA-seq variants are true variants.
			recordStatus = self._from_RNA(record)
		elif training_type == 'VAF':
			recordStatus = self._from_VAF(record)
		else:
			message = "The training type is not supported: '{}'".format(self.truthset_type)
			raise ValueError(message)

		return recordStatus

	def _combineTruthsets(self, samples):
		""" Combines the truthsets of several samples into a single file.
			Parameters
			----------
				samples: list<dict<>>
					A list of outputs from self._per_sample()
						* 'filename-indel'
						* 'filename-snv'
			Example Picard Command
			java -jar picard.jar SortVcf \
				I=vcf_1.vcf \
				I=vcf_2.vcf \
				O=sorted.vcf
		"""
		truthset_folder = getPipelineFolder('truthset', 'multiple')
		sampleIds = sorted(i['PatientID'] for i in samples)
		# Extracts patient id from the barcode. ex. TCGA-2H-A9GF -> A9GF
		sampleIds = [i.split('-')[-1] for i in sampleIds]
		sampleIds = ",".join(sampleIds)

		indel_truthset_filename = os.path.join(
			truthset_folder, "{}.{}.merged.indel.truthset.vcf".format(
				sampleIds,
				training_type
			)
		)
		snp_truthset_filename = os.path.join(
			truthset_folder, "{}.{}.merged.snp.truthset.vcf".format(
				sampleIds,
				training_type
			)
		)

		# Generate a command-line to parse the files.
		# snp_cmd_string   = "--variant " + " --variant ".join([i['filename-snv']   for i in samples]) + '\\'
		# indel_cmd_string = "--variant " + " --variant ".join([i['filename-indel'] for i in samples]) + '\\'
		snp_cmd_string   = " ".join(["I={}".format(i['filename-snv'])   for i in samples])
		indel_cmd_string = " ".join(["I={}".format(i['filename-indel']) for i in samples])
		
		picard_base_command = """java -jar {program} SortVcf {variants} O={output}"""
		snv_command = picard_base_command.format(
			program     = self.picard_program,
			reference   = self.reference,
			variants    = indel_cmd_string,
			output      = indel_truthset_filename)

		indel_command = picard_base_command.format(
			program     = self.picard_program,
			reference   = self.reference,
			variants    = snp_cmd_string,
			output      = snp_truthset_filename)

		systemtools.Terminal(snv_command,   show_output = True)
		systemtools.Terminal(indel_command, show_output = True)

		result = {
			'PatientID': sampleIds,
			'filename-indel': indel_truthset_filename,
			'filename-snv': snp_truthset_filename
		}

		return result

	#@staticmethod
	def _from_intersection(self, record, **kwargs):

		_separator = '-'

		# GATK includes callers that filtered the variant, designated by a "filterIn{caller}" label.
		callset = [i for i in record.INFO['set'].split(_separator) if 'filter' not in i.lower()]

		# MuSE and Somaticsniper are snp-only, and should not be counted in the intersection of indel callsets
		if record.is_indel:
			num_callers_in_intersection = self.indel_intersection
		else:
			num_callers_in_intersection = self.snp_intersection

		_is_intersection = len(callset) == 1 and callset[0] == 'Intersection'

		_is_in_n_sets = len(callset) >= num_callers_in_intersection
		
		validation_status = _is_intersection or _is_in_n_sets

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
			'validation status': 0 if filterOut is False else 1
		}
		return row

	def _from_VAF(self, record):
		sample_vaf = self._getSampleVAF(record)
		normal_vaf = sample_vaf['NORMAL']['vaf']
		tumor_vaf = sample_vaf['TUMOR']['vaf']

		# Filter out variants that were rejected by a caller
		# This is only needed for variants in the callsets of a single caller
		filter_out = False
		# Validate variants according to VAF status.
		# This is only for testing purposes, a better version should be used later.
		high_tumor_vaf  = tumor_vaf  >= self.min_tumor_vaf
		high_normal_vaf = normal_vaf >= self.max_normal_vaf
		_somatic_vaf = (high_tumor_vaf and not high_normal_vaf)
		_somatic_vaf = _somatic_vaf or (not high_tumor_vaf and normal_vaf == 0.0)
		if _somatic_vaf and not filter_out:
			validation_status = 1  # Somatic
		elif (not high_tumor_vaf and (normal_vaf > 0.0 and not high_normal_vaf)) or filter_out:
			validation_status = 0  # Non-Somatic
		else: validation_status = 2  # 'Unknown'

		validation_status = int(validation_status == 1)

		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'VAF',
			'validation status': validation_status
		}
		return row
	
	def export(self):
		result = {
			'indel': self.indel_filename,
			'snp':   self.snp_filename
		}
		return result


class VCFtoMAF:
	def __init__(self, sample, input_vcf, options):
		""" Converts a VCF file to and annotated MAF file. Include additional fields such as VAF.
			Parameters
			----------
				sample:
				input_vcf: string, dict<string:string>
					The vcf file to convert.
				options: dict<>
			Additional Columns
			------------------
				VAF: The variant allele frequency

		"""
		self.reference_build    = 'GRCh38'
		self.vcftomaf_program   = options['Programs']['vcf2maf']
		self.vep_program        = options['Programs']['varianteffectpredictor']
		self.reference          = options['Reference Files']['reference genome']
		self.maf_fieldnames     = [ 
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
			"CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "ALLELE_NUM"
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
		patientId = sample['PatientID']
		output_folder = getPipelineFolder('vcftomaf', patientId)
		raw_maf = os.path.join(output_folder, "{0}.raw.maf".format(patientId))
		self.maf = os.path.join(output_folder, "{0}.maf".format(patientId))
		maf_file = self.vcftomaf(sample, input_vcf, raw_maf)
		print(maf_file)
		print("Output MAF: ", self.maf)

	@staticmethod
	def modifyMAF(sample, input_maf, truthset):
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

				row['Tumor_Sample_Barcode']         = sample['SampleID']
				row['Matched_Norm_Sample_Barcode']  = sample['NormalID']
				row['tumor_bam_uuid']               = sample['SampleUUID']
				row['normal_bam_uuid']              = sample['NormalUUID']
				validation_status = truthset(
					sample = sample['SampleID'],
					chrom = row['Chromosome'],
					pos = row['Start_Position'])

				row['Validation_Status_VAF']        = validation_status.get('VAF')
				row['Validation_Status_RNA']        = validation_status.get('RNA-seq')
				row['Validation_Status_Intersect']  = validation_status.get('Intersection')

				row['Filter_Status'] = 0

				output_maf.append(row)

		return output_maf

	def vcftomaf(self, sample, input_vcf, output_maf):
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
				script      = self.vcftomaf_program,
				vep         = self.vep_program,
				reference   = self.reference,
				ncbi        = self.reference_build,
				tumor       = sample['SampleID'],
				normal      = sample['NormalID'],
				vcf         = input_vcf,
				maf         = output_maf)
		if not os.path.isfile(output_maf):
			systemtools.Terminal(command)

		return output_maf

	@staticmethod
	def filterMAF(sample, options, input_maf, caller):
		"""
			Filter Criteria
			---------------
				dbSNP: filter out variants present in dbSNP
				Indel Filtering: filter out variants in known indels
				Germline/normal filtering: filter out variants present in a panel of normals. (derive from GDC MAFs?)
				PON Filtering
				
		"""
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
	def __init__(self, training_samples, prediction_samples, options, **kwargs):
		""" Parameters
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
			print("Pipeline: training samples: ")
			for element in training_samples:
				print("\t", element['PatientID'])
			print()
			print("Pipeline: prediction samples: ")
			for element in prediction_samples:
				print("\t", element['PatientID'])
		
		##################### Pre-process the callsets of the training samples. ####################
		for sample in training_samples + prediction_samples:
			#### Fix the files that will be used to generate the truthset (the training samples) ###
			print("Pipeline: Processing ", sample['PatientID'])
			print("\tPipeline: retrieving the original callset")
			original_callset     = getSampleCallset(sample['PatientID'], 'original')
			fixed_callset_folder = getPipelineFolder('callset', sample['PatientID'], 'original-fixed')

			somaticseq_folder = options['Programs']['SomaticSeq']
			print("\tvcftools.fixCallerOutputs()")
			fixed_callset = vcftools.fixCallerOutputs(
				original_callset, somaticseq_folder, output_folder = fixed_callset_folder)
			################### Separate the corrected files into indels and snps ##################
			print("\tgetCallsetFolder()")
			split_callset_output_folder = getPipelineFolder(
				'callset', sample['PatientID'], 'original-fixed-split')
			print("\tvcftools.splitCallset()")
			split_callset = vcftools.splitCallset(fixed_callset, split_callset_output_folder)

		# Generate a truthset using the training samples. Separate indel/snp truthsets will be generated.
		training_sample = training_samples[0]

		snp_intersection_requirement = 4
		indel_intersection_requirement = 2
		
		truthset    = Truthset(
			training_sample, 
			options, 
			truthset_type      = 'Intersection', 
			snp_intersection   = snp_intersection_requirement,
			indel_intersection = indel_intersection_requirement,
			**kwargs)

		#somaticseq  = SomaticSeqPipeline(training_samples, prediction_samples, options, truthset = truthset)
		# Build the model using the truthset and training samples.
		somaticseq_trained_model = SomaticSeqPipeline(
			kind = 'trainer',
			sample = training_sample, 
			options = options,
			truthset = truthset.export()['snp']
			)
		snp_classifier_filename = somaticseq_trained_model.export()['classifier']
		# Run the prediction model.
		somaticseq_prediction_model = SomaticSeqPipeline(
			kind = 'prediction',
			sample = prediction_samples[0],
			options = options,
			truthset = None,
			snp_classifier = snp_classifier_filename)

		# list of dictionaries with the output filenames for a given sample.
		#predictions = somaticseq.predictions

		# Filter
		"""
		filtered_samples = list()
		for sample_vcfs in predictions:
			sample_vcf = [i for i in prediction_samples if i['PatientID'] == sample_vcfs['PatientID']][0]
			indel_vcf = sample_vcf['indelVCF']
			snv_vcf   = sample_vcf['snpVCF']
			sample_vcf['indelFilteredVCF'] = FilterVCF(indel_vcf, options).output
			sample_vcf['snvFilteredVCF']   = FilterVCF(snv_vcf, options).output

			filtered_samples.append(sample_vcfs)
		
		# Annotate and convert to MAF
		maf_samples = list()
		for filtered_sample in filtered_samples:
			sample = [i for i in prediction_samples if i['PatientID'] == filtered_sample['PatientID']]
			indel_vcf = filtered_sample['indelFilteredVCF']
			snv_vcf = filtered_sample['snvFilteredVCF']
			indel_maf = VCFtoMAF(sample, indel_vcf, options)
			snv_maf = VCFtoMAF(sample, snv_vcf, options)
			filtered_sample['indelMAF'] = indel_maf
			filtered_sample['snvMAF'] = snv_maf
			maf_samples.append(filtered_sample)
		"""
		# Combine MAFs


"""
	Truthset -> Somaticseq training -> somaticseq prediction -> filter -> convert to MAF - > Combine MAFs
"""
###################################################################################################
########################## Define global objects used in the pipeline #############################
###################################################################################################
options = configparser.ConfigParser()
options.read(OPTIONS_FILENAME)

#GET_CALLSET = callertools.CallerOutputClassifier(reduce = True, verbose = True)
GET_CALLSET = callertools.CallerClassifier()
GATK_MERGE_CALLSET = GATKMergeSampleCallsets(options)

###################################################################################################
###################################### Run the Pipeline ##########################################
###################################################################################################


if __name__ == "__main__" and True:
	documents_folder = os.path.join(os.getenv('HOME'), 'Documents')

	full_sample_list_filename = os.path.join(documents_folder, "DNA-Seq_sample_list.tsv")

	full_sample_list = tabletools.readCSV(full_sample_list_filename)

	training_type = 'Intersection'

	training_sample_ids = ['TCGA-2H-A9GR']
	prediction_sample_ids = ['TCGA-L5-A43M']#, 'TCGA-R6-A6Y2-CHR1']

	training_samples    = [i for i in full_sample_list if i['PatientID'] in training_sample_ids]
	prediction_samples  = [i for i in full_sample_list if i['PatientID'] in prediction_sample_ids]

	Pipeline(training_samples, prediction_samples, options)

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
	print(sample)
elif True:
	pass
