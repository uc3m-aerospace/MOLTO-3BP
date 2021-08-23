import sys, os, yaml
sys.path.insert(0, os.getcwd())

from Common.Base.singleton import SingletonMeta


class Configurator(metaclass=SingletonMeta):
    """
    Class
    """
    config = {}
    working_path = ""
    current_config_file = ""

    def __init__(self):
        self.working_path = os.getcwd() + "/Config"
        self.config_file = self.working_path + "/prod.configuration.yaml"
        self.read_config()

    def get_config_file(self):
        """
        Load the configuration
        """
        if os.path.isfile(self.working_path + "/local.configuration.yaml"):
            self.config_file = self.working_path + "/local.configuration.yaml"
        elif os.path.isfile(self.working_path + "/dev.configuration.yaml"):
            self.config_file = self.working_path + "/dev.configuration.yaml"

        if not os.path.isfile(self.config_file):
            raise FileNotFoundError('Config file not found')
        return self.config_file

    def read_config(self):
        """
        Parse the configuration YAML format
        """
        with open(self.get_config_file()) as stream:
            self.config = yaml.safe_load(stream)

    def get(self, variable, default=None):
        """
        Get variable from configuration
        """
        parts = variable.split('.')
        value = None
        config = self.config
        for key in parts:
            if key in config:
                value = config[key]
                config = config[key]
            else:
                value = None
        return value if value is not None else default
