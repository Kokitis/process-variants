import os

DOCUMENTS_FOLDER = os.path.join(os.getenv('HOME'), 'Documents')

if os.name == 'nt':
	GITHUB_FOLDER = os.path.join(os.getenv('USERPROFILE'), 'Documents', 'Github')
else:
	GITHUB_FOLDER = os.path.join(DOCUMENTS_FOLDER, 'Github')

PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
OPTIONS_FILENAME = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_configuration.txt")

options = configparser.ConfigParser()
options.read(OPTIONS_FILENAME)

full_sample_list_filename = os.path.join(documents_folder, "DNA-Seq_sample_list.tsv")

full_sample_list = tabletools.readCSV(full_sample_list_filename)