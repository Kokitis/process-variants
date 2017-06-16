
import os
import configparser
import global_variables as gv

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
		self.vcftomaf_program   = gv.options['Programs']['vcf2maf']
		self.vep_program        = gv.options['Programs']['varianteffectpredictor']
		self.reference          = gv.options['Reference Files']['reference genome']
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
		patientId 		= sample['PatientID']
		output_folder 	= getPipelineFolder('vcftomaf', patientId)
		raw_maf 		= os.path.join(output_folder, "{0}.raw.maf".format(patientId))
		self.maf 		= os.path.join(output_folder, "{0}.maf".format(    patientId))
		maf_file 		= self.vcftomaf(sample, input_vcf, raw_maf)
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


if __name__ == "__main__":
	pass