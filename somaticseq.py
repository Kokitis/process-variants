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
import pytools.filetools as filetools
import pytools.tabletools as tabletools
import pytools.timetools as timetools

import varianttools.vcftools as vcftools
import varianttools.callertools as callertools


def getBasename(filename):
	bn = os.path.basename(filename)
	bn = os.path.splitext(bn)[0]
	return bn

class SomaticSeqPipeline:
	def __init__(self, sample, output_folder, options):
		self.output_folder = output_folder
		self.somaticseq_folder = options['Programs']['SomaticSeq']
		self.somaticseq_program = os.path.join(self.somaticseq_folder, "SomaticSeq.Wrapper.sh")
		self.modify_vjsd_script = os.path.join(self.somaticseq_folder, "modify_VJSD.py")
		self.ada_trainer_script = os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_builder.R")
		self.ada_prediction_script=os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_predictor.R")
		self.tsv_to_vcf_script  = os.path.join(self.somaticseq_folder, "SSeq_tsv2vcf.py")
		self.merged_vcf2tsv_script     = os.path.join(self.somaticseq_folder, "SSeq_merged.vcf2tsv.py")
		self.gatk_program = options['Programs']['GATK']

		training_folder = self._getSomaticseqFolder('trainer', 'TCGA-2H-A9GR')
		self.snp_classifier = os.path.join(
			training_folder, 
			"TCGA-2H-A9GR.training.modified.merged.excluded.snp.tsv.Classifier.RData")

		self.targets    = sample['ExomeTargets']
		self.cosmic     = options['Reference Files']['cosmic']
		self.dbSNP      = options['Reference Files']['dbSNP']
		self.reference  = options['Reference Files']['reference genome']

		self.truthset = os.path.join(
			PIPELINE_FOLDER, 
			"truthset", 
			"files", 
			"TCGA-2H-A9GR-CHR1", 
			"TCGA-2H-A9GR-CHR1.Intersection.snp.final.truthset.vcf"
		)

		self.runWorkflow(sample, callset)

	def _convertToTable(self, sample, callset, merged_callset, output_folder):
		print("Converting to a table...")
		timer = timetools.Timer()
		output_filename = os.path.join(
			output_folder,
			"{}.training.modified.merged.excluded.snp.tsv".format(sample['PatientID'])
		)
		print("Merged_callset: {}\t{}".format(os.path.exists(merged_callset), merged_callset))
		command = """python3 {script} \
			-ref {reference} \
			-nbam {normal} \
			-tbam {tumor} \
			-truth {truth} \
			-dbsnp {dbSNP} \
			-cosmic {cosmic} \
			-myvcf {merged} \
			-muse {muse} \
			-mutect {mutect2} \
			-sniper {somaticsniper} \
			-strelka {strelka} \
			-varscan {varscan} \
			-outfile {output}""".format(
				script 		= self.merged_vcf2tsv_script,
				reference 	= self.reference,
				normal = sample['NormalBAM'],
				tumor  = sample['TumorBAM'],
				dbSNP  = self.dbSNP,
				cosmic = self.cosmic,
				truth  = self.truthset,
				merged = merged_callset,
				muse   			= callset['muse-snp'],
				mutect2 		= callset['mutect2-snp'],
				somaticsniper 	= callset['somaticsniper-snp'],
				strelka 		= callset['strelka-snp'],
				varscan 		= callset['varscan-snp'],
				output 			= output_filename
			)
		#if not os.path.exists(output_filename):
		systemtools.Terminal(command, use_system = True)
		print("\tResult: {}\t{}".format(os.path.exists(output_filename), output_filename))
		print(timer.to_iso())
		return output_filename

	def _concatPaths(self, *paths):
		if isinstance(paths[0], (tuple, list)):
			paths = paths[0]
		paths = [i for i in paths if i]
		path = os.path.join(self.output_folder, *paths)
		return path

	def _getSomaticseqFolder(self, step, patientId):
		""" Retrieves the path to a folder reserved for a specific
			step in the pipeline.
			Pipeline Structure:
				callsets
				truthset
					files
						
				somaticseq
					prediction
						[training_type]
					training
						[training_type]
				vcftomaf
		"""
		if step == 'files':
			folder_names = 'files'
		elif step == 'callset':
			folder_names = ()
		elif step == 'somaticseq':
			folder_names = ('somaticseq')
		elif step == 'trainer':
			folder_names = ('training', patientId)
		elif step == 'prediction':
			folder_names = ('prediction', patientId)
		else:
			folder_names = None
		if folder_names is None:
			message = "The indicated pipeline step is not defined: " + step
			raise KeyError(message)
		else:
			pipeline_folder = self._concatPaths(folder_names)

		filetools.checkDir(pipeline_folder, True)
		return pipeline_folder

	def _mergeVariantFiles(self, callset, output_folder, patientId):
		print("Merging files...")
		output_filename = os.path.join(
			output_folder,
			"{}.training.modified.merged.vcf".format(patientId)
		)

		command = """java -jar {gatk} \
			-T CombineVariants \
			-R {reference} \
			-genotypeMergeOptions UNSORTED \
			-V {muse} \
			-V {mutect2} \
			-V {somaticsniper} \
			-V {strelka} \
			-V {varscan} \
			-o {output}""".format(
			gatk = self.gatk_program,
			reference = self.reference,
			muse = 		callset['muse-snp'],
			mutect2 = 	callset['mutect2-snp'],
			somaticsniper = callset['somaticsniper-snp'],
			strelka = 	callset['strelka-snp'],
			varscan = 	callset['varscan-snp'],
			output = 	output_filename
		)
		if not os.path.exists(output_filename):
			systemtools.Terminal(command, use_system = True)
		print("\tResult: {}\t{}".format(os.path.exists(output_filename), output_filename))
		return output_filename

	def runWorkflow(self, sample, callset):
		print("Running Workflow...")
		patientId = sample['PatientID']
		method = 'prediction'
		classifier = callertools.CallerClassifier()

		original_callset_folder = "/home/upmc/Documents/Genomic_Analysis/1_callsets/{}/original_callset".format(
			sample['PatientID'])

		split_callset_folder = "/home/upmc/Documents/Genomic_Analysis/1_callsets/{}/original_corrected_split_callset".format(
			sample['PatientID'])

		original_callset = classifier(original_callset_folder)
		vcftools.splitCallset(original_callset, split_callset_folder)
		callset = classifier(split_callset_folder, type = 'snp')
		pprint(callset)
		process_output_folder = self._getSomaticseqFolder(method, patientId)
		
		processed_callset = self._preprocessVariantFiles(
			callset, 
			process_output_folder, 
			patientId)
		
		merged_raw_variant_file = self._mergeVariantFiles(
			processed_callset, 
			process_output_folder, 
			patientId)
		
		reduced_raw_variant_file= self._reduceVariantTargets(
			merged_raw_variant_file, 
			process_output_folder, 
			patientId)
		
		trained_snp_table = self._convertToTable(
			sample, 
			callset, 
			reduced_raw_variant_file, 
			process_output_folder)

		if method == 'trainer':
			self.buildTrainer(trained_snp_table)
		elif method == 'prediction':
			self.runPredictor(trained_snp_table)

		return trained_snp_table

	def buildTrainer(self, input_filename):
		print("Building model...")
		timer = timetools.Timer()
		command = "{script} {infile}".format(
			script = self.ada_trainer_script,
			infile = input_filename
		)
		systemtools.Terminal(command, use_system = True)
		print(timer.to_iso())

	def runPredictor(self, input_filename):

		"""
		"""
		basename = getBasename(input_filename)
		output_folder = os.path.dirname(input_filename)
		output_filename = "{}.predicted_scores.tsv".format(basename)
		output_filename = os.path.join(output_folder, output_filename)

		print("Predicting Scores...")
		timer = timetools.Timer()
		command = "{script} {classifier} {infile} {outfile}".format(
			script = self.ada_prediction_script,
			classifier = self.snp_classifier,
			infile = input_filename,
			outfile = output_filename
		)
		systemtools.Terminal(command, use_system = True)
		print(timer.to_iso())

	def _reduceVariantTargets(self, input_filename, output_folder, patientId):
		print("Excluding non-exome targets...")
		output_filename = os.path.join(
			output_folder,
			"{}.training.modified.merged.excluded.vcf".format(patientId)
		)

		command = """intersectBed -header -a {infile} -b {targets} > {output}""".format(
			infile = input_filename,
			targets = self.targets,
			output = output_filename
		)

		systemtools.Terminal(command, use_system = True)
		print("\tResult: {}\t{}".format(os.path.exists(output_filename), output_filename))
		return output_filename


	def _preprocessVariantFiles(self, callset, output_folder, patientId):
		print("Processing Variant Files...")
		processed_callset = dict()
		for caller_name, filename in callset.items():
			processed_file = self._processVJSDFiles(caller_name, filename, output_folder, patientId)
			processed_callset[caller_name] = processed_file

		return processed_callset

	def _processVJSDFiles(self, caller, input_file, output_folder, patientId):
		"""
			Parameters
			----------
				input_file:
				output_folder:
				caller: {'varscan', 'somaticsniper', 'muse'}
		"""
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
			command = "python3 {program} -method {method} -infile {infile} -outfile {outfile}".format(
				program = self.modify_vjsd_script,
				method = method,
				infile = input_file,
				outfile = output_filename)
			systemtools.Terminal(command, use_system = True)
		else:
			shutil.copy2(input_file, output_filename)			

		return output_filename

PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")
options = configparser.ConfigParser()
options.read(OPTIONS_FILENAME)

if __name__ == "__main__":
	#sample_id = "TCGA-2H-A9GR"
	sample_id = "TCGA-L5-A43M"

	full_sample_list_filename = os.path.join(
		os.path.dirname(PIPELINE_FOLDER),
		"DNA-Seq_sample_list.tsv")

	full_sample_list = tabletools.Table(full_sample_list_filename)

	sample = full_sample_list('PatientID', sample_id)

	_callset_folder = "/home/upmc/Documents/Genomic_Analysis/1_callsets/TCGA-2H-A9GR/original_callset/"


	callset = {
		'muse': 		_callset_folder + 'TCGA-L5-A43M-11A_vs_TCGA-L5-A43M-01A.Muse.chr1.vcf',
		'mutect2': 		_callset_folder + 'TCGA-2H-A9GR-11A_vs_TCGA-2H-A9GR-01A.mutect2.chr1.vcf',
		'somaticsniper':_callset_folder + 'TCGA-2H-A9GR-11A_vs_TCGA-2H-A9GR-01A.somaticsniper.chr1.vcf',
		'strelka-snp': 	_callset_folder + 'TCGA-2H-A9GR-11A_vs_TCGA-2H-A9GR-01A.passed.somatic.snvs.vcf.strelka.chr1.vcf',
		'varscan-snp':	_callset_folder + 'TCGA-2H-A9GR-11A_vs_TCGA-2H-A9GR-01A.raw.snp.Somatic.hc.chr1.vcf'
	}
	test_output_folder = "/home/upmc/Documents/Genomic_Analysis/debug_folder"
	SomaticSeqPipeline(sample, test_output_folder, options)


