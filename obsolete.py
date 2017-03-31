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
