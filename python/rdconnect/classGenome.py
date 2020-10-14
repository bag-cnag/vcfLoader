"""Classes

This module contains the abstract definition of class 'GenomicData' and the two 
functions used in 'vcfLoader' to attach new methods to it.

This file should be imported as module to create new 'GenomicData' objects 
and/or to attach new methods to the class.

To add new methods to the class:
	* 'add_method' allows to attach a single method and give it a name
	* 'add_funcs_from_module' allows to attach all public functions of a module
"""

from inspect import getmembers, isfunction

class GenomicData():
	"""
	A class used to represent any type of genomic data through the pipeline
	implemented in 'vcfLoader'.

	It is an abstract class with no attributes and no methods, only '__str__' is
	overloaded. Methods and attributes are being attached to the class 
	dynamically using the function 'add_method' and 'add_funcs_from_module'.
	"""
	def __str__(self):
		return '<GenomicData>'


def add_method(function, calling):
	"""Adds a new method to GenomicData

	This functions creates a new method in 'GenomicData' class with the name 
	'calling' using the implementation in 'function'.

	Parameters
	----------
	function: func, mandatory
		The function to be attached to 'GenomicData' class.
	calling: str, mandatory
		Name of the new created method.

	Raises
	------
	Exception
		If 'function' is not a function or 'calling' is not a string.
	"""
	if type(calling) is not str:
		raise Exception('Given a "calling" that is no string.')
	if not callable(function):
		raise Exception('Given a "function" that can not be called.')
	setattr(GenomicData, calling, function)


def add_funcs_from_module(module):
	"""Adds all the public functions of a module as new method to GenomicData

	This functions creates a new method in 'GenomicData' class for each of the 
	public function in the given module. A function is considered public if it
	does not start with underscore ('_').

	Parameters
	----------
	module: module, mandatory
		The module from witch to get all functions to be added to 'GenomicData'
	
	Raises
	------
	Exception
		If the modules has items that can be understood as function but are
		not functions such as variables.
	"""
	fun_list = [x for x in getmembers(module) if isfunction(x[1])]
	pub_fun = [x for x in fun_list if not x[0].startswith('_')]
	for nam, fun in pub_fun:
		add_method(nam, fun)

