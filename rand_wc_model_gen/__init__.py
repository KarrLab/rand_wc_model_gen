import pkg_resources

# read version
with open(pkg_resources.resource_filename('rand_wc_model_gen', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

# API
from . import config
from . import kb_gen
from . import model_gen
