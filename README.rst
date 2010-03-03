About drosophiline
=================

This will be a collection of scripts
which do application-specific
pre-processing and post-processing
around hidden Markov state inference.


Requirements
============

Operating system requirements:

* This project was developed using Ubuntu,
  so it will probably work on Debian-based Linux distributions.
* It might work with non-Debian-based Unix variants.
* It probably will not work on Windows.

Major dependencies:

* A recent version of Python-2.x_ (2.6+).
* A C compiler which is not too different from gcc.

Python package and module dependencies:

* numpy_
* argparse_


Installation
============

First install ``virtualenv`` and ``pip``.
Next activate a virtual environment.

Drosophiline can be installed directly from its github_
repository as follows::

    $ pip install git+git://github.com/argriffing/drosophiline

Alternatively if you have cloned the git repo
as ``~/repos/drosophiline`` for some other reason,
drosophiline can be installed from this local repository as follows::

    $ pip install -e ~/repos/drosophiline


Uninstallation
--------------

Uninstall drosophiline using::

    $ pip uninstall drosophiline


Usage
=====

For now, the project will consist only of
scripts run from the command line.


.. _Python-2.x: http://www.python.org
.. _argparse: http://code.google.com/p/argparse
.. _virtualenv: http://virtualenv.openplans.org
.. _pip: http://pip.openplans.org
.. _pypi: http://pypi.python.org
.. _github: http://github.com
.. _numpy: http://numpy.scipy.org
