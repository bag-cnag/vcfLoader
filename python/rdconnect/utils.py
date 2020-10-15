import os

def chrom_str_to_int(chrom):
	"""This function is used parse a chromosome from 'string name' to 'numeric
	name'.

	This functions converts string names for chromosome to numeric names 
	according to HAIL numerical system where:

		* X: 24
		* Y: 25
		* MT: 23

	Parameters
	----------
	chrom: str, mandatory
		A chromosome identifier ('1', '2', 'X', 'MT', ...).
	
	Returns
	-------
	A string where the string names is changed to numeric name.
	"""
	if chrom in ('X', 'Y', 'MT'):
		return {'X': '24', 'Y': '25', 'MT': '23'}[chrom]
	if chrom.lower() == 'all':
		return ''
	return chrom


def create_chrom_filename(filename, chrom):
	"""This function is used to create names using CNAG-CRG's BU nomenclature.

	This functions is sued to create filenames using CNAG-CRG's BU nomenclature
	where the templates contains the tag 'chromosome' and it has to be replaced
	by the chromosome's numeric name.

	Parameters
	----------
	filename: str, mandatory
		Template used to create the chromosome's iterative files.
	chrom: str, mandatory
		A chromosome's numeric name.
	
	Returns
	-------
	The original name for the file with the right chromosome in it.
	"""
	return filename.replace('chromosome', str(chrom))


def destination_dense_matrices(destination_path, filename, somatic = False):
	"""This function returns the path to a file annotated with VEP.

	Parameters
	----------
	destination_path: str, mandatory
		Path where the file will be saved.
	filename: str, mandatory
		Template used to create the chromosome's iterative files.
	somatic: bool, mandatory
		Indicates if the saved wile contains somatic mutations (set it to True)
		or germline (set it to False).

	Returns
	-------
	A string with the path to save the annotated VEP file.
	"""
	return os.path.join(destination_path, 'dense_matrices{}'.format('_somatic' if somatic else ''), filename)


def destination_vep(destination_path, filename, somatic = False):
	"""This function returns the path to a file annotated with VEP.

	Parameters
	----------
	destination_path: str, mandatory
		Path where the file will be saved.
	filename: str, mandatory
		Template used to create the chromosome's iterative files.
	somatic: bool, mandatory
		Indicates if the saved wile contains somatic mutations (set it to True)
		or germline (set it to False).

	Returns
	-------
	A string with the path to save the annotated VEP file.
	"""
	return os.path.join(destination_path, 'annotated_vep{}'.format('_somatic' if somatic else ''), filename)

def destination_dbnsfp(destination_path, filename, somatic = False):
	"""This function returns the path to a file annotated with VEP.

	Parameters
	----------
	destination_path: str, mandatory
		Path where the file will be saved.
	filename: str, mandatory
		Template used to create the chromosome's iterative files.
	somatic: bool, mandatory
		Indicates if the saved wile contains somatic mutations (set it to True)
		or germline (set it to False).

	Returns
	-------
	A string with the path to save the annotated VEP file.
	"""
	return os.path.join(destination_path, 'annotated_dbnsfp{}'.format('_somatic' if somatic else ''), filename)