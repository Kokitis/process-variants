import os
import sys
import configparser
import shutil
from pprint import pprint

if os.name == 'nt':
	GITHUB_FOLDER = os.path.join(os.getenv('USERPROFILE'), 'Documents', 'Github')
else:
	GITHUB_FOLDER = os.path.join(os.getenv('HOME'), 'Documents', 'Github')
sys.path.append(GITHUB_FOLDER)

import pytools.systemtools as systemtools

import pytools.tabletools as tabletools
import pytools.timetools as timetools

import varianttools.vcftools as vcftools
import varianttools.callertools as callertools
import time
from pipeline_manager import getPipelineFolder

def getBasename(filename):
	bn = os.path.basename(filename)
	bn = os.path.splitext(bn)[0]
	return bn

class SomaticSeqPipeline:
	def __init__(self, kind, sample, options, truthset = None, snp_classifier = None, options = None):
		"""
			Parameters
			----------
				kind: {'trainer', 'prediction', 'table'}
					* 'trainer': Trains the classifier.
						Requires 'truthset'
					* 'prediction': Calculated the probability that a variant is a true variant.
						Requires 'snp_classifier', the output from the 'trainer' step.
					* 'table': generates the table used for the training and prediction steps.

		"""
		if options is None: options = default_options
		print("Runing Somaticseq...")
		print("\tMode: ", kind)
		print("\tsample: ", sample['PatientID'])

		self.kind = kind
		self.somaticseq_folder 		= options['Programs']['SomaticSeq']
		self.somaticseq_program 	= os.path.join(self.somaticseq_folder, "SomaticSeq.Wrapper.sh")
		self.modify_vjsd_script 	= os.path.join(self.somaticseq_folder, "modify_VJSD.py")
		self.ada_trainer_script 	= os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_builder.R")
		self.ada_prediction_script	= os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_predictor.R")
		self.tsv_to_vcf_script  	= os.path.join(self.somaticseq_folder, "SSeq_tsv2vcf.py")
		self.merged_vcf2tsv_script  = os.path.join(self.somaticseq_folder, "SSeq_merged.vcf2tsv.py")
		self.gatk_program = options['Programs']['GATK']

		self.targets    = sample['ExomeTargets']
		self.cosmic     = options['Reference Files']['cosmic']
		self.dbSNP      = options['Reference Files']['dbSNP']
		self.reference  = options['Reference Files']['reference genome']

		self.snp_classifier = snp_classifier # Should only be provided in training mode.

		self.truthset = truthset # Should be None if in prediction/table mode

		self.runWorkflow(sample)
	
	def export(self):
		result = {
			'classifier': self.classifier,
			'table': self.trained_snp_table
		}
		return result
	
	def runWorkflow(self, sample):
		print("Running Workflow...")
		patientId = sample['PatientID']
		callset = self._getRawCallset(patientId)
		pprint(callset)
		process_output_folder = getPipelineFolder('somaticseq-' + self.kind, patientId)
		
		processed_callset = self._processVJSDFiles(
			callset, 
			process_output_folder, 
			patientId
		)
		
		merged_raw_variant_file = self._mergeVariantFiles(
			processed_callset, 
			process_output_folder, 
			patientId
		)
		
		reduced_raw_variant_file= self._reduceVariantTargets(
			merged_raw_variant_file, 
			process_output_folder, 
			patientId
		)
		
		self.trained_snp_table = self._convertToTable(
			sample, 
			callset, 
			reduced_raw_variant_file, 
			process_output_folder
		)

		if self.kind == 'trainer':
			self.classifier = self.buildTrainer(self.trained_snp_table)
		elif self.kind == 'prediction':
			self.classifier = None
			prediction_table = self.runPredictor(self.trained_snp_table)
			prediction_vcf = self._convertToVcf(sample, prediction_table)
		elif kind == 'table':
			pass

		return self.trained_snp_table

	def _getRawCallset(self, patientId):
		classifier = callertools.CallerClassifier()

		original_callset_folder = getPipelineFolder('callset', patientId, 'original')

		split_callset_folder = getPipelineFolder('callset', patientId, 'original-split')

		original_callset = classifier(original_callset_folder)
		
		vcftools.splitCallset(original_callset, split_callset_folder)
		
		callset = classifier(split_callset_folder, type = 'snp')	

		return callset

	def _processVJSDFiles(self, callset, output_folder, patientId):
		"""
			Parameters
			----------
				input_file:
				output_folder:
				caller: {'varscan', 'somaticsniper', 'muse'}
		"""
		processed_callset = dict()
		for caller, input_file in callset.items():
			print("\tProcessing {}...".format(caller))
			if caller == 'muse':
				method = 'MuSE'
			elif caller == 'somaticsniper':
				method = 'SomaticSniper'
			elif 'varscan' in caller:
				method = 'VarScan2'
			else:
				method = None

			output_filename = os.path.join(
				output_folder, 
				"{}.training.modified.vcf".format(patientId)
			)

			if method:
				command = """python3 {program} \
					--call-method {method} \
					--input-vcf {infile} \
					--output-vcf {outfile}""".format(
					program = self.modify_vjsd_script,
					method = method,
					infile = input_file,
					outfile = output_filename)
				systemtools.Terminal(command, use_system = True)
			else:
				shutil.copy2(input_file, output_filename)

			processed_callset[caller] = output_filename


		return processed_callset

	def _mergeVariantFiles(self, callset, output_folder, patientId):
		"""
			Note: The only records that are merged are those that are
			unfiltered in at least one caller.
		"""
		print("Merging files...")
		output_filename = os.path.join(
			output_folder,
			"{}.{}.modified.merged.vcf".format(patientId, self.kind)
		)

		command = """java -jar {gatk} \
			--analysis_type CombineVariants \
			--reference_sequence {reference} \
			--genotypemergeoption UNSORTED \
			--filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
			--variant {muse} \
			--variant {mutect2} \
			--variant {somaticsniper} \
			--variant {strelka} \
			--variant {varscan} \
			--out {output}""".format(
			gatk 		= self.gatk_program,
			reference 	= self.reference,
			muse 		= callset['muse-snp'],
			mutect2 	= callset['mutect2-snp'],
			somaticsniper = callset['somaticsniper-snp'],
			strelka 	= callset['strelka-snp'],
			varscan 	= callset['varscan-snp'],
			output 		= output_filename
		)

		if not os.path.exists(output_filename):
			systemtools.Terminal(command, use_system = True)
		return output_filename
	
	def _reduceVariantTargets(self, input_filename, output_folder, patientId):
		print("Excluding non-exome targets...")
		output_filename = os.path.join(
			output_folder,
			"{}.{}.modified.merged.excluded.vcf".format(patientId, self.kind)
		)

		command = """intersectBed -header -a {infile} -b {targets} > {output}""".format(
			infile = input_filename,
			targets = self.targets,
			output = output_filename
		)

		systemtools.Terminal(command, use_system = True)
		#print("\tResult: {}\t{}".format(os.path.exists(output_filename), output_filename))
		return output_filename
	
	def _convertToTable(self, sample, callset, merged_callset, output_folder):
		print("Converting to a TSV file...")
		start_time = time.time()
		output_filename = os.path.join(
			output_folder,
			"{}.{}.modified.merged.excluded.snp.tsv".format(sample['PatientID'], self.kind)
		)
		print("Merged_callset: {}\t{}".format(os.path.exists(merged_callset), merged_callset))
		command = """python3 {script} \
			--p-scale phred \
			--genome-reference {reference} \
			--normal-bam-file {normal} \
			--tumor-bam-file {tumor} \
			--dbsnp-vcf {dbSNP} \
			--cosmic-vcf {cosmic} \
			--vcf-format {merged} \
			--muse-vcf {muse} \
			--mutect-vcf {mutect2} \
			--somaticsniper-vcf {somaticsniper} \
			--strelka-strelka-vcf {strelka} \
			--varscan-vcf {varscan} \
			--output-tsv-file {output}""".format(
				script 			= self.merged_vcf2tsv_script,
				reference 		= self.reference,
				normal 			= sample['NormalBAM'],
				tumor  			= sample['TumorBAM'],
				dbSNP  			= self.dbSNP,
				cosmic 			= self.cosmic,
				truth  			= self.truthset,
				merged 			= merged_callset,
				muse   			= callset['muse-snp'],
				mutect2 		= callset['mutect2-snp'],
				somaticsniper 	= callset['somaticsniper-snp'],
				strelka 		= callset['strelka-snp'],
				varscan 		= callset['varscan-snp'],
				output 			= output_filename
			)

		if self.kind == 'trainer':
			command += " --ground-truth-vcf " + self.truthset

		#if not os.path.exists(output_filename):
		if not os.path.exists(output_filename):
			systemtools.Terminal(command, use_system = True)
		else:
			print("The Somaticseq table already exists.")
		print("\tResult: {}\t{}".format(os.path.exists(output_filename), output_filename))
		stop_time = time.time()

		duration = timetools.Duration(seconds = stop_time - start_time)
		print("Converted the table in ", duration.isoformat())
		return output_filename

	def _convertToVcf(self, sample, input_filename):
		print("Converting to a VCF file...")

		output_filename = os.path.splitext(input_filename)[0] + ".vcf"

		command = """python3 {script} \
			--tsv-in {infile} \
			--vcf-out {outfile} \
			--normal-sample-name {normalid} \
			--tumor-sample-name {tumorid} \
			--emit-all \
			--individual-mutation-tools {tools} \
			--phred-scale""".format(
				script = self.tsv_to_vcf_script,
				infile = input_filename,
				outfile = output_filename,
				normalid = sample['NormalID'],
				tumorid = sample['SampleID'],
				tools = "MuSE CGA SomaticSniper Strelka VarScan2"
			)
		print(command)
		if not os.path.exists(output_filename):
			systemtools.Terminal(command, use_system = True)
		return output_filename
		

	def buildTrainer(self, input_filename):
		print("Building model...")
		command = "{script} {infile}".format(
			script = self.ada_trainer_script,
			infile = input_filename
		)
		expected_output = input_filename + '.Classifier.RData'
		print("expected_output: ", expected_output)
		print(command)
		if not os.path.exists(expected_output):
			systemtools.Terminal(command, use_system = True)
		else:
			print("The Somaticseq classifier already exists.")
		return expected_output

	def runPredictor(self, input_filename):

		"""
		"""
		basename = getBasename(input_filename)
		output_folder = os.path.dirname(input_filename)
		output_filename = "{}.predicted_scores.tsv".format(basename)
		output_filename = os.path.join(output_folder, output_filename)

		print("Predicting Scores...")
		command = "{script} {classifier} {infile} {outfile}".format(
			script = self.ada_prediction_script,
			classifier = self.snp_classifier,
			infile = input_filename,
			outfile = output_filename
		)
		print(command)
		systemtools.Terminal(command, use_system = True)

		return output_filename



PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")
default_options = configparser.ConfigParser()
default_options.read(OPTIONS_FILENAME)

if __name__ == "__main__":
	#sample_id = "TCGA-2H-A9GR"
	sample_id = "TCGA-L5-A43M"

	full_sample_list_filename = os.path.join(
		os.path.dirname(PIPELINE_FOLDER),
		"DNA-Seq_sample_list.tsv")

	full_sample_list = tabletools.Table(full_sample_list_filename)

	sample = full_sample_list('PatientID', sample_id)

	test_output_folder = "/home/upmc/Documents/Genomic_Analysis/debug_folder"
	SomaticSeqPipeline(sample, test_output_folder, default_options)



