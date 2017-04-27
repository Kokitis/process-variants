import os
import collections
import sys
import configparser
import math

import vcf
import csv
from pprint import pprint

"""
	1. Harmonize VCF fields
	2. Merge caller VCFs
	3. Convert VCFs to MAF
	4. Combine Patient MAFs
"""
PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
CONSOLE_LOG_FILE = os.path.join(PIPELINE_FOLDER, "system_log_file.txt")
OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")

TEST_MAF_FILE = ""
TEST_SEG_FILE = ""

def Terminal(command, useSystem = True, toFile = False):
	if useSystem:
		os.system(command)
		output = "os.system({0})".format("command")
	else:
		result = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		output = str(result.stdout.read(),'utf-8')

	if toFile:
		with open(CONSOLE_LOG_FILE, 'a') as console_log_file:
			console_log_file.write(output)

	return output

def checkdir(path, full = False):
	if not os.path.isdir(path):
		if full:	os.makedirs(path)
		else:		os.mkdir(path)

def readTSV(filename, headers = False):
	with open(filename, 'r') as file1:
		reader = csv.DictReader(file1, delimiter = '\t')
		fieldnames = reader.fieldnames
		reader = list(reader)

	if headers: return reader, fieldnames
	else:       return reader

def writeTSV(table, filename, fieldnames = None):
	if fieldnames is None:
		fieldnames = sorted(table[0].keys())
	with open(filename, 'w', newline = '') as file1:
		writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = fieldnames)
		writer.writeheader()
		writer.writerows(table)
	return filename

class Workflow:
	def __init__(self, **kwargs):
		self.setupEnvironment(**kwargs)
		self.runWorkflow(**kwargs)
	def runCommand(self, command, expected_output):

		if not os.path.exists(expected_output):
			Terminal(command)
		else:
			pass

		return os.path.exists(expected_output)

class Gistic:
	"""
		Requirements
		------------
			segmentation file
				Contains the segmented data identified by CBS.
				(1)  Sample           	(sample name)
				(2)  Chromosome  		(chromosome number)
				(3)  Start Position  	(segment start position, in bases)
				(4)  End Position  		(segment end position, in bases)
				(5)  Num markers      	(number of markers in segment)
				(6)  Seg.CN       		(log2() -1 of copy number)
			reference genome
				The reference genome.
		Optional
		--------
			markers file
				Identifies marker positions used for array or capture experiements.
				The easiest way to create a markers file, is to take the first 3 columns of the CN file used as input to the segmentation algorithm that produced your seg file.
		cnv file
			identified germilne cnvs to exclude from the analysis.
	"""
	def __init__(self, options, seg_file, refgene, markers = None):
		print("Running GISTIC with unvetted options!")
		self.gistic_folder = options['Programs']['GISTIC']
		self.gistic_program = os.path.join(self.gistic_folder, 'gistic2')
		self.refgene = os.path.join(self.gistic_folder, 'refgenefiles', 'hg38.UCSC.add_miR.160920.refgene.mat')
		self.output_folder = ""
		self.reference = options['Reference Files']['reference']

		#seg_file = self._convertSegmentationFile(seg_file)
		
		self.output = self.run(seg_file, markers)

	def _convertSegmentationFile(self, seg_file):
		""" Converts the log ratios of the segmentation file into log2-1 (such that a copynumber of 2 is 0) """
		output_file = os.path.splitext(seg_file)[0] + '.log2.tsv'
		segments, fieldnames = readTSV(seg_file, True)
		table = list()
		for row in segments:
			mean = float(row['Segment_Mean'])
			if mean == 0.0:
				row['Segment_Mean'] = -5.0
			else:
				row['Segment_Mean'] = math.log2(mean) - 1
			table.append(row)
		writeTSV(table, output_file, fieldnames)

		return output_file


	def run(self, seg_file, markers = None):

		"""./gistic2 -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -alf $alf -cnv $cnvfile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme"""
		if markers is None: markers = 10
		command = """{program} \
			-b {output_folder} \ 
			-seg {segfile} \
			-mk {markers} \ 
			-refgene {reference} 
			-genegistic 1 \
			-smallmem 1 \
			-broad 1 \
			-brlen 0.5 \ 
			-conf 0.90 \
			-armpeel 1 \
			-savegene 1 \
			-gcm extreme""".format(
				program = self.gistic_program,
				output_folder = self.output_folder,
				reference = self.reference,
				markers = markers,
				segfile = seg_file
			)
		Terminal(command)
		expected_output = {}
		return expected_output

	def test(self):
		seg_file = TEST_SEG_FILE
		refgene = "/home/upmc/Documents/Reference/Genome/GRCh38.d1.vd1.fa"
		markers = None

		g = cls(seg_file, refgene, markers)

		return g


class MutSigCV(Workflow):
	def setupEnvironment(self, kwargs):
		reference_files_folder = ""
		self.output_folder = ""
		self.matlab_program = ""
		self.coverage_file = os.path.join(reference_files_folder, "exome_full192.coverage.txt")
		self.covariate_file = os.path.join(reference_files_folder, "gene.covariates.txt")

	def runWorkflow(self, **kwargs):
		"""
			Keyword Arguments
			-----------------
				'maf_file'
		"""
		#command = "run_MutSigCV.sh <path_to_MCR> mutations.maf coverage.txt covariates.txt output.txt"
		basename = os.path.basename(maf_file)
		output_file = os.path.join(self.output_folder, basename + '.mutsig.sig_genes.txt')
		command = "run_MutSigCV.sh {matlab} {maf_file} {coverage} {covariates} {output}".format(
			matlab 		= self.matlab_program,
			maf_file 	= maf_file,
			coverage 	= self.coverage_file,
			covariates 	= self.covariates_file,
			output_file = output_file)

		status = self.runCommand(command, output_file)
		return status

class SciClone:
	""" 
		Requires
		--------
			seg_files (need to convert log-ratio to absolute values with 2*2^log-ratio?)
			vaf_files: Tab delimited file detaiiling the vaf at each position.
				Format: chr, pos, ref_reads, var_reads, vaf
	"""
	def __init__(self, seg_files, vaf_files):
		pass

	def generateRScript(self, seg_files, vaf_files):

		r_script = """library(sciClone)
			#read in vaf data from three related tumors
			#format is 5 column, tab delimited: 
			#chr, pos, ref_reads, var_reads, vaf

			v1 = read.table("data/vafs.tumor1.dat",header=T);
			v2 = read.table("data/vafs.tumor2.dat",header=T);
			v3 = read.table("data/vafs.tumor3.dat",header=T);

			#read in regions to exclude (commonly LOH)
			#format is 3-col bed
			regions = read.table("data/exclude.loh")

			#read in segmented copy number data
			#4 columns - chr, start, stop, segment_mean   
			cn1 = read.table("data/copy_number_tum1")
			cn2 = read.table("data/copy_number_tum2")
			cn3 = read.table("data/copy_number_tum3")

			#set sample names
			names = c("Sample1","Sample2","Sample3")


			#Examples:
			#------------------------------------
			#1d clustering on just one sample
			sc = sciClone(vafs=v1,
			         copyNumberCalls=cn1,
			         sampleNames=names[1],
			         regionsToExclude=reg1)
			#create output
			writeClusterTable(sc, "results/clusters1")
			sc.plot1d(sc,"results/clusters1.1d.pdf")

			#------------------------------------
			#2d clustering using two samples:
			sc = sciClone(vafs=list(v1,v2),
			              copyNumberCalls=list(cn1,cn2),
			              sampleNames=names[1:2],
			               regionsToExclude=regions)
			#create output
			writeClusterTable(sc, "results/clusters2")
			sc.plot1d(sc,"results/clusters2.1d.pdf")
			sc.plot2d(sc,"results/clusters2.2d.pdf")


			#------------------------------------
			#3d clustering using three samples:
			sc = sciClone(vafs=list(v1,v2,v3),
			              copyNumberCalls=list(cn1,cn2,cn3),
			              sampleNames=names[1:3],
			               regionsToExclude=regions)
			#create output
			writeClusterTable(sc, "results/clusters2")
			sc.plot1d(sc,"results/clusters2.1d.pdf")
			sc.plot2d(sc,"results/clusters2.2d.pdf")
			sc.plot3d(sc, sc@sampleNames, size=700, outputFile="results/clusters3.3d.gif")

			#This pattern generalizes up to N samples, except for plotting, which caps out at 3d for obvious reasons."""


			r_script = '\n'.join(r_script.splitlines())

def getCopynumberVariants(folder, caller = ""):
	""" Retrieves all copynumber files for the project.
	"""
	file_list = list()
	if caller == 'Varscan':
		_format = "varscan.copynumber.called.segments"
	elif caller == 'CNVkit':
		_format = ""
	else:
		_format = ""

	for barcode in os.listdir(folder):
		caller_path = os.path.join(folder, barcode, caller)
		segment_file = search(caller_path, _format)

	return file_list



if __name__ == "__main__" and True:
	options = configparser.ConfigParser()
	options.read(OPTIONS_FILENAME)
	
	project_maf = ""
	project_segments = ""

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