[//]: # ( [![PyPI package](https://img.shields.io/pypi/v/rand_wc_model_gen.svg)](https://pypi.python.org/pypi/rand_wc_model_gen) )
[![Documentation](https://readthedocs.org/projects/rand_wc_model_gen/badge/?version=latest)](http://docs.karrlab.org/rand_wc_model_gen)
[![Test results](https://circleci.com/gh/KarrLab/rand_wc_model_gen.svg?style=shield)](https://circleci.com/gh/KarrLab/rand_wc_model_gen)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/rand_wc_model_gen/badge.svg)](https://coveralls.io/github/KarrLab/rand_wc_model_gen)
[![Code analysis](https://api.codeclimate.com/v1/badges/a9d32ece26a8d3c363e0/maintainability)](https://codeclimate.com/github/KarrLab/rand_wc_model_gen)
[![License](https://img.shields.io/github/license/KarrLab/rand_wc_model_gen.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/rand_wc_model_gen/README.md?pixel)

# rand_wc_model_gen

The whole-cell model generator generates synthetic whole-cell model descriptions that are used to test downstream components of the whole-cell (WC) modeling pipeline. In particular, these synthetic models can be used to test model description systems and model simulators.

To test modeling pipeline components, synthetic WC models are better than WC models of real cells for several reasons:

* Synthetic models are scalable
* Synthetic models can have unusual properties
* Synthetic models can be generated quickly

This initial version of the model generator creates reaction network models.

## Installation
1. Install dependencies
2. Install this package 
  ```
  pip install git+git://github.com/KarrLab/rand_wc_model_gen#egg=rand_wc_model_gen
  ```

## Documentation
Please see the [API documentation](http://docs.karrlab.org/rand_wc_model_gen).

## License
The build utilities are released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Karr Lab](http://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

## Questions and comments
Please contact the [Karr Lab](http://www.karrlab.org) with any questions or comments.