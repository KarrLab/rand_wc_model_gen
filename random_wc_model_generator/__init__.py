import pkg_resources

# read version
with open(pkg_resources.resource_filename('random_wc_model_generator', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()