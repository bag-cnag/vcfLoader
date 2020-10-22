"""classConfig

This module contains the implementation of the class 'ConfigFile' and its
auxiliary functions for its creation.

This file should be imported as module to create new 'ConfigFile' object, 
containing all the configuration for any step in the annotation pipeline.

How to create a ConfigFile?

	>>> x = ConfigFile('/path/to/config/json/config.json')

How to get a configuration parameter?

	>>> x['section/parameter']
	>>> x['section/subsection/parameter']

Do you need to overwrite a parameter to match the expected behaviour of the
pipeline?

	>>> x = x.overwrite('section/parameter', 'new_value')

"""

# from classConfig import ConfigFile
# x = ConfigFile('/Users/chernandez/Desktop/test_config.json')

import json
import operator
from copy import deepcopy
from functools import reduce
from traceback import format_exc


def _get_from_dict(data, key_list):
	return reduce(operator.getitem, key_list, data)


def _get_keys_dict(dic, parent = None):
	rst = []
	for k, v in dic.items():
		if isinstance(v, dict):
			if parent is None:
				rst += _get_keys_dict(v, k)
			else:
				rst += _get_keys_dict(v, parent + '/' + k)
		else:
			if parent is None:
				rst.append(k)
			else:
				rst.append(parent + '/' + k)
	return rst


def _set_nested(dic, key_list, value):
	for key in key_list[:-1]:
		dic = dic.setdefault(key, {})
	dic[key_list[-1]] = value


class ConfigFile:

	def __init__(self):
		self.config_path = None
		self.data = {}

	def __getitem__(self, key):
		if key in self.keys:
			return _get_from_dict(self.data, key.split('/'))
		else:
			return None

	def __setitem__(self, key, value):
		_set_nested(self.data, key.split('/'), value)
		self.keys = _get_keys_dict(self.data)

	def __str__(self):
		return 'ConfigFile<{}>'.format(str(self.config_path))

	def get(self, key):
		return self[key]

	def overwrite(self, key, value):
		new = deepcopy(self)
		_set_nested(new.data, key.split('/'), value)
		new.keys = _get_keys_dict(new.data)
		return new

	@staticmethod
	def from_file(config_path):
		x = ConfigFile()
		x.config_path = config_path
		x.data = {}
		with open(config_path) as data_file:
			x.data = json.load(data_file)
		try:
			x.keys = _get_keys_dict(x.data)
		except Exception as ex:
			raise Exception('Invalid configuration file provided.\n{}\n{}'.format(str(ex), str(format_exc())))
		return x

	@staticmethod
	def from_str(data_str):
		x = ConfigFile()
		x.config_path = None
		x.data = json.loads(data_str)
		try:
			x.keys = _get_keys_dict(x.data)
		except Exception as ex:
			raise Exception('Invalid configuration file provided.\n{}\n{}'.format(str(ex), str(format_exc())))
		return x
