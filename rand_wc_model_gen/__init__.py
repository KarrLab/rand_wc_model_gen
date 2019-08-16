import pkg_resources

# read version
with open(pkg_resources.resource_filename('rand_wc_model_gen', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

# API
# from . import analysis
# from . import config
# from wc_kb_gen import random

from .model_gen.core import RandModelGen
