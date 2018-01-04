[//]: # ( [![PyPI package](https://img.shields.io/pypi/v/wc_model_generator.svg)](https://pypi.python.org/pypi/wc_model_generator) )
[![Documentation](https://readthedocs.org/projects/wc_model_generator/badge/?version=latest)](http://wc_model_generator.readthedocs.org)
[![Test results](https://circleci.com/gh/KarrLab/wc_model_generator.svg?style=shield)](https://circleci.com/gh/KarrLab/wc_model_generator)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/wc_model_generator/badge.svg)](https://coveralls.io/github/KarrLab/wc_model_generator)
[![Code analysis](https://codeclimate.com/github/KarrLab/wc_model_generator/badges/gpa.svg)](https://codeclimate.com/github/KarrLab/wc_model_generator)
[![License](https://img.shields.io/github/license/KarrLab/wc_model_generator.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/wc_model_generator/README.md?pixel)

# wc_model_generator

The whole-cell model generator generates synthetic whole-cell model descriptions that are used to test downstream components of the whole-cell (WC) modeling pipeline. In particular, these synthetic models can be used to test model description systems and model simulators.

Synthetic WC models are better for testing modeling pipeline components than WC models of real cells for several reasons:

* Synthetic models are scalable
* Synthetic models can have unusual properties
* Synthetic models can be generated quickly

This initial version of the model generator creates reaction network models.

## Installation
1. Install dependencies
2. Install this package 
  ```
  pip install git+git://github.com/KarrLab/wc_model_generator#egg=wc_model_generator
  ```

## Documentation
Please see the [API documentation](http://wc_model_generator.readthedocs.io).

## License
The build utilities are released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Karr Lab](http://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

## Questions and comments
Please contact the [Karr Lab](http://www.karrlab.org) with any questions or comments.