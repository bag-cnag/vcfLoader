import json
import operator
from copy import deepcopy
from functools import reduce

# from classConfig import ConfigFile
# x = ConfigFile('/Users/chernandez/Desktop/test_config.json')
#

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

	def __init__(self, config_path):
		self.config_path = config_path
		self.data = {}
		with open(config_path) as data_file:
			self.data = json.load(data_file)
		self.keys = _get_keys_dict(self.data)

	def __getitem__(self, key):
		if key in self.keys:
			return _get_from_dict(self.data, key.split('/'))
		else:
			return None

	def __str__(self):
		return 'ConfigFile<{}>'.format(self.config_path)

	def get(self, key):
		return self[key]

	def overwrite(self, key, value):
		new = copy.deepcopy(self)
		_set_nested(new.data, key.split('/'), value)
		return new
