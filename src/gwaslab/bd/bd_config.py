from os import path
from pathlib import Path

class Options_dic:
    def __init__(self,path_dic):
        self.paths = path_dic
        self.default = path_dic.copy()
    def set_option(self,key,path):
        self.paths[key] = path
    def reset_option(self):
        self.paths = self.default.copy()

#default_dic={
#    "config": path.dirname(__file__) + "/data/config.json",
#    "reference": path.dirname(__file__) + "/data/reference.json",
#    "formatbook": path.dirname(__file__) + "/data/formatbook.json",
#    "data_directory": path.expanduser('~') + "/.gwaslab/"
#}

default_dic={
    "config": path.join(Path(__file__).parents[1], "data","config.json"),
    "reference": path.join(Path(__file__).parents[1], "data","reference.json"),
    "formatbook": path.join(Path(__file__).parents[1], "data","formatbook.json"),
    "data_directory": path.expanduser('~') + "/.gwaslab/"
}


options = Options_dic(default_dic)