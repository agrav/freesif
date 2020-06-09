# FREESIF

Work easily with data from SESAM Interface Files.

## General

### About

Python package that enables easy extraction of data from SESAM Interface Files.

The following file formats are currently supported:

- Formatted finite element model file (.FEM)
- Formatted interface file (.SIF)
- Unformatted interface file (.SIU)

### Getting started

Install the latest release

```console
pip install freesif
```

... and import it into a script

```python
from freesif import open_sif
```

### Resources

* [**Source**](https://github.com/agrav/freesif)
* [**Issues**](https://github.com/agrav/freesif/issues)
* [**Changelog**](https://github.com/agrav/freesif/releases)
* [**Documentation**](https://github.com/agrav/freesif/blob/master/README.md)
* [**Download**](https://pypi.org/project/freesif/)

## Contribute

These instructions will get you a copy of the project up and running on your local machine for development and testing
purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Install Python version 3.6 or later from either https://www.python.org or https://www.anaconda.com.

### Clone the source code repository

At the desired location, run:

```git clone https://github.com/agrav/freesif.git```

### Installing

To get the development environment running:

... create an isolated Python environment and activate it,

```console
python -m venv /path/to/new/virtual/environment

/path/to/new/virtual/environment/Scripts/activate
```

... install the dev dependencies in [requirements.txt](requirements.txt),

```console
pip install -r requirements.txt
```

.. and install the package in development mode.

```console
python setup.py develop
```

You should now be able to import the package in the Python console,

```python
import freesif
help(freesif)
```

### Running the tests

The unit tests are automated using the `unittest` and `pytest` framework. Run the test by...

```console
pytest --cov=freesif --cov-report term-missing tests/
```

### Building the package

Build tarball and wheel distributions by:

```console
pip install wheel
python setup.py sdist bdist_wheel
```

The distribution file names adhere to the [PEP 0427](https://www.python.org/dev/peps/pep-0427/#file-name-convention)
convention `{distribution}-{version}(-{build tag})?-{python tag}-{abi tag}-{platform tag}.whl`.

<!---
### Building the documentation

The html documentation is build using [Sphinx](http://www.sphinx-doc.org/en/master)

```console
sphinx-build -b html docs\source docs\_build
```
--->

### Deployment
Packaging, unit testing and deployment to [PyPi](https://pypi.org) is automated using GitHub Actions.

### Versioning

We apply the "major.minor.micro" versioning scheme defined in [PEP 440](https://www.python.org/dev/peps/pep-0440/).

Cut a new version by applying a Git tag like `1.0.1` at the desired commit and then
[setuptools_scm](https://github.com/pypa/setuptools_scm/#setup-py-usage) takes care of the rest. For the versions
available, see the [tags on this repository](https://github.com/agrav/freesif/tags).

## Authors

* **Audun Gravdal Johansen** - [agrav](https://github.com/agrav)

## Maintainers

* **Audun Gravdal Johansen** - [agrav](https://github.com/agrav)
* **Per Voie** - [tovop](https://github.com/tovop)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

