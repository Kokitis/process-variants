
import os
import pandas
import re
import shutil
import vcf
from pprint import pprint

def readCSV(filename, headers = False, **kwargs):
	delimiter = kwargs.get('delimiter', kwargs.get('sep', '\t'))
	with open(filename, 'r') as file1:
		reader = csv.DictReader(file1, delimiter = delimiter)
		fieldnames = reader.fieldnames
		reader = list(reader)

	if headers: return reader, fieldnames
	else:       return reader

def writeCSV(table, filename, fieldnames = None):
	if fieldnames is None:
		fieldnames = sorted(table[0].keys())
	with open(filename, 'w', newline = '') as file1:
		writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = fieldnames)
		writer.writeheader()
		writer.writerows(table)
	return filename

def checkDir(path, force = False):
	if os.path.isfile(path):
		folder, basename = os.path.split(path)
	else:
		folder = path

	if not os.path.exists(folder):
		if force:
			os.makedirs(folder)
		else:
			os.mkdir(folder)
	return os.path.exists(folder)

PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
OPTIONS_FILENAME = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")

def compareCallers(left, right):
	""" Uses GATK CombineVariants to merge and compare the output of two callers.
 		Parameters
 		----------
 			left: string
 			right: string
 	"""
	gatk_program = options['Programs']['GATK']
	reference = options['Reference Files']['reference genome']
	output_file = ""

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
		-priority {rod}""".format(
			gatk = gatk_program,
			reference = reference,
			left = left,
			right = right,
			output = output_file,
			rod = order)
def countVariants(filename):
	""" Counts the number of variants detected per Chromosome."""
	chromosomes = dict()
	with open(filename, 'r') as file1:
		reader = vcf.Reader(file1)

		for record in reader:
			if record.CHROM not in chromosomes:
				chromosomes[record.CHROM] = 1
			else:
				chromosomes[record.CHROM] += 1

	pprint(chromosomes)

def combineDepthOfCoverageOBS(filenames):
	regex = re.compile("[-A-Z]+")
	barcodes = dict()
	gene_list = None
	for filename in filenames:
		path, basename = os.path.split(filename)
		barcode = regex.search(basename)

		patient_table = pandas.read_csv(filename, sep = '\t')
		if gene_list is None:
			gene_list = patient_table['Genes'].values
		barcodes[barcode] = patient_table['average_coverage'].values

	coverage_table = [gene_list]
	fieldnames = ['Genes']
	for barcode, coverage in barcodes.items():
		fieldnames.append(barcode)
		coverage_table.append(coverage)

def combineDepthOfCoverage(filenames):
	tables = list()
	for filename in filenames:
		data = pandas.read_csv(filename, delimiter = "\t")
		data.set_index('Genes', inplace = True)
		tables.append(data['average_coverage'])

	table = pandas.DataFrame(tables)
	return table

	writeCSV(coverage_table, "", fieldnames)



class SearchOutputFiles:
	def __init__(self, source, destination = None, DEBUG = False):
		self.DEBUG = DEBUG

		self.callers = self._defineFilenameTemplates()

		matches = dict()
		for barcode in os.listdir(source):
			barcode_folder = os.path.join(source, barcode)
			matches[barcode] = self.search(barcode_folder)
		self.matches = matches
		if destination is not None:
			self.copyCaseFolders(destination, matches)
		#pprint(matches)

	def __call__(self, barcode):

		return self.matches.get(barcode, [])

	@staticmethod
	def _compile(regex):
		sample_regex = "[-a-z0-9]+"
		prefix_regex = "{0}_vs_{0}".format(sample_regex)
		regex = {k: re.compile(v.format(sample = sample_regex, prefix = prefix_regex)) for k, v in regex.items()}
		return regex

	def _getNewFilename(key, basename):
		""" Renames the output files. """
		basename, ext = os.path.splitext(basename)
		if key in {'STR0', 'STR1', 'STR2', 'STR3'}:
			filename = basename + '.strelka.vcf'
		elif key in {}:
			pass

		return filename

	def _defineFilenameTemplates(self):


		cnvkit = {
			'CNV0': "{sample}.cnr",
			'CNV1': "{sample}.cns",
			'CNV2': "{sample}.targetcoverage",
			'CNV3': "{sample}.antitargetcoverage"
		}

		depthofcoverage = {
			'DOC0': "{prefix}.depthofcoverage"
		}

		freec = {
			'FRE0': "{prefix}.freec"
		}

		haplotypecaller = {
			'HAP0': "{prefix}.haplotypecaller"
		}

		muse = {
			'MUS0': "{prefix}.muse.vcf"
		}

		mutect2 = {
			'MUT0': "{prefix}.mutect2"
		}

		somaticsniper = {
			'SOM0': "{prefix}.somaticsniper.lq",
			'SOM1': "{prefix}.somaticsniper.hq",
			'SOM2': "{prefix}.somaticsniper.vcf"
		}

		strelka = {
			'STR0': "all.somatic.indels",
			'STR1': "all.somatic.snvs",
			'STR2': "passed.somatic.indels",
			'STR3': "passed.somatic.snvs"
		}

		unifiedgenotyper = {
			'UNI0': "{prefix}.unifiedgenotyper"
		}

		varscan = {
			'VAR0': "{prefix}.varscan"
		}

		cnvkit 			= self._compile(cnvkit)
		depthofcoverage = self._compile(depthofcoverage)
		freec 			= self._compile(freec)
		haplotypecaller = self._compile(haplotypecaller)
		muse 			= self._compile(muse)
		mutect2 		= self._compile(mutect2)
		somaticsniper 	= self._compile(somaticsniper)
		strelka 		= self._compile(strelka)
		unifiedgenotyper= self._compile(unifiedgenotyper)
		varscan 		= self._compile(varscan)

		original_names = {
			'cnvkit': cnvkit,
			'depthofcoverage': depthofcoverage,
			'freec': freec,
			'haplotypecaller': haplotypecaller,
			'muse': muse,
			'mutect2': mutect2,
			'somaticsniper': somaticsniper,
			'strelka': strelka,
			'unifiedgenotyper': unifiedgenotyper,
			'varscan': varscan
		}
		return callers

	def copyCaseFolders(self, destination, matches):
		for barcode, filenames in matches.items():
			destination_folder = os.path.join(destination, barcode)
			if len(filenames) == 0: continue
			checkDir(destination_folder)
			for filename in filenames:
				path, basename = os.path.split(filename)
				destination_file = os.path.join(destination_folder, basename.lower())
				shutil.copyfile(filename, destination_file)

	def isOutput(self, abs_path):
		path, basename = os.path.split(abs_path)

		for caller, regexes in self.callers.items():
			for key, regex in regexes.items():
				#print(key, regex)
				match = regex.search(basename.lower())
				if self.DEBUG:
					print(key, '\t', match)
				if match is not None:
					return caller
		else:
			return None

	def search(self, source):
		""" """
		matches = list()
		for item in os.listdir(source):
			abs_path = os.path.join(source, item)
			if os.path.isdir(abs_path):
				matches += self.search(abs_path)
			else:
				match = self.isOutput(abs_path)
				if match is not None:
					matches.append(abs_path)
		return matches


if __name__ == "__main__":
	filename = "C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Data\\Genomic_Analysis\\1_input_vcfs\\"
	destination = "C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Data\\Genomic_Analysis\\1_input_vcfs\\test\\"
	SearchOutputFiles(filename, destination)