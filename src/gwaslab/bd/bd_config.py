from os import path
from pathlib import Path

class Options_dic:
    """
    A class to manage and modify configuration paths for gwaslab.

    Attributes:
        paths (dict): Current configuration paths
        default (dict): Default paths for reset functionality
    """
    
    def __init__(self, path_dic):
        """
        Initialize with a dictionary of configuration paths.

        Parameters:
            path_dic (dict): Dictionary containing initial path configurations
        """
        self.paths = path_dic
        self.default = path_dic.copy()
    
    def set_option(self, key, path):
        """
        Update a specific configuration path.

        Parameters:
            key (str): Configuration key to update
            path (str): New path value for the specified key
        """
        self.paths[key] = path
    
    def reset_option(self):
        """
        Reset all paths to their default values.
        """
        self.paths = self.default.copy()

    def get_option(self):
        """
        Return the options as a dict
        """
        return self.paths

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
