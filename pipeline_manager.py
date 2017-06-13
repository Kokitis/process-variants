"""
	Manages the file structure of the pipeline.
"""
import os
import filetools.filetools as filetools
import varianttools.callertools as callertools
PIPELINE_FOLDER = "/home/upmc/Documents/Genomic_Analysis"
GET_CALLSET = callertools.CallerClassifier()
def _concatPaths(*paths):
	if isinstance(paths[0], (tuple, list)):
		paths = paths[0]
	paths = [i for i in paths if i]
	path = os.path.join(PIPELINE_FOLDER, *paths)
	return path

def getCallsetFolder(patientId, kind):
	if kind == 'original':
		subfolder = "original_callset"
	elif kind == 'original-split':
		subfolder = 'original_split_callset'
	elif kind == 'original-fixed':
		subfolder = 'original_corrected_callset'
	elif kind == 'original-fixed-split':
		subfolder = "original_corrected_split_callset"
	elif kind == 'merged':
		subfolder = 'merged_callset'
	else:
		subfolder = None

	if subfolder is not None:
		subfolder = ['1_callsets', patientId, subfolder]

	return subfolder


def getPipelineFolder(step, patientId, kind = None, **kwargs):
	step = step.split('-')
	if 'callset' in step:
		subfolders = getCallsetFolder(patientId, kind)
	elif 'truthset' in step:
		step.remove('truthset')
		subfolders = _getTruthsetFolder(patientId)
	elif 'somaticseq' in step:
		step.remove('somaticseq')
		kind = step[0]
		subfolders = _getSomaticseqFolder(patientId, kind)
	else:
		subfolders = None

	if subfolders is None:
		message = "'{}' is not a valid step".format(step)
		raise ValueError(message)

	folder = _concatPaths(*subfolders)
	filetools.checkDir(folder, True)
	return folder

def _getSomaticseqFolder(patientId, kind):
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
	if kind == 'trainer':
		folder_names = ['somaticseq', 'trainer', patientId]
	elif kind == 'prediction':
		folder_names = ['somaticseq', 'prediction', patientId]
	else:
		folder_names = None

	return folder_names

def _getTruthsetFolder(patientId):
	#subfolders = [training_type, patientId]
	#subfolders = ['truthset'] + subfolders
	subfolders = ['truthset', patientId]
	return subfolders


def getSampleCallset(patientid, callset_type, group = 'all'):
	""" Search and select a specific callset for a sample.
		Parameters
		----------
			patientid:  string, dict<>
				The patient barcode or a dictionary containing the patient barcode saved
				under the key 'PatientID'.
			callset_type: {'original', 'original-fixed', 'original-fixed-split'}
			group: {'all', 'indel', 'snp'}; default 'all
		Returns
		-------
			A dictionary mapping callers to output files.
	"""
	if isinstance(patientid, dict):
		patientid = patientid['PatientID']

	# Terms to determine which folders to skip when searching for variant files.
	_exclusion_terms = ['strelka', 'chromosome']

	callset_folder = getCallsetFolder(patientid, kind = callset_type)
	sample_variants = GET_CALLSET(
		callset_folder, 
		kind = callset_type, 
		exclude = _exclusion_terms, 
		logic = 'and'
	)

	if group == 'all':
		pass
	else:
		if group not in {'indel', 'snp', 'snv'}:
			message = "The provided callset group does not exist: " + group
			raise ValueError(message)
		elif 'split' not in callset_type:
			message = "The callset is not split into indels or snps: " + callset_type
			raise ValueError(message)
		sample_variants = {k.split('-')[0]: v for k, v in sample_variants.items() if group in k}

	return sample_variants


