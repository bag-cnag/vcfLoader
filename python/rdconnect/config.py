import json

def readConfig(config_path):
    with open(config_path) as data_file:
        data = json.load(data_file)
    return data

def readFilesList(list_path):
    with open(list_path) as files_paths:
        paths = files_paths.read().splitlines()
        print paths
    return paths
