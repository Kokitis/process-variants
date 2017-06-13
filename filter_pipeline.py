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
		}

		return result

	@staticmethod
	def writeVCF(records, filename, reader):
		with open(filename, 'w') as file1:
			writer = vcf.Writer(file1, reader)
			for record in records:
				writer.write_record(record)
		return filename

	def _getFilterStatus(self, record):
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

		filter_status = dbSNP_status or known_indel_status
		filter_status = filter_status or base_quality_status
		return filter_status

	def _filterFromdbSNP(self, record):
		message = "_filterFromdbSNP is not implemented."
		raise NotImplementedError(message)

	def _filterKnownIndels(self, record):
		message = "_filterKnownIndels is not implemented."
		raise NotImplementedError(message)

	def _filterBaseQuality(self, record):
		message = "_filterBasequality is not implemented."
		raise NotImplementedError(message)

	def _filterConsequence(self, record):
		message = "_filterConsequence is not implemented."
		raise NotImplementedError(message)

	def _filterGermlineVariants(self, record):
		message = "_filterGermlineVariants is not implemented."
		raise NotImplementedError(message)
