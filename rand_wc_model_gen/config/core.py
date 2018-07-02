""" Configuration

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import configobj
import os
import pkg_resources
import rand_wc_model_gen
import wc_utils.config.core


def get_config(extra_path=None, extra_vals=None):
    """ Get configuration

    Args:
        extra_path (:obj:`dict`, optional): additional configuration to override
        extra_vals (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.config.core.ConfigPaths(
        default=pkg_resources.resource_filename('rand_wc_model_gen', 'config/core.default.cfg'),
        schema=pkg_resources.resource_filename('rand_wc_model_gen', 'config/core.schema.cfg'),
        user=[
            'rand_wc_model_gen_without_git.cfg',
            os.path.expanduser('~/.wc/rand_wc_model_gen.cfg'),
        ],
    )

    if extra_path:
        paths.user.insert(0, extra_path)

    context = {
        'package_path': pkg_resources.resource_filename('rand_wc_model_gen', ''),
        'version': rand_wc_model_gen.__version__,
    }
    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra_vals, context=context)
