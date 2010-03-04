"""
This is the setup script.
This script is automatically run by
easy_install or by pip on the user's machine when
she installs the module from pypi.

Here is some documentation for this process.
http://docs.python.org/extending/building.html

More info:
http://wiki.python.org/moin/Distutils/Tutorial

Register the metadata with pypi as follows:
python setup.py register

Send to pypi as follows:
python setup.py sdist upload --show-response
"""

from distutils.core import setup
from distutils.core import Extension

from drosophiline import drosophiline_version


scripts = [
        'bin/dline-create-transition-params',
        'bin/dline-create-emission-params',
        'bin/dline-compute-likelihoods',
        'bin/dline-call-polymorphisms',
        'bin/dline-user-friendly']

classifiers = [
        'Development Status :: 1 - Planning',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics']
        

# This seems like a standard method.
long_description = open('README.rst').read()

setup(
        name = 'drosophiline',
        version = drosophiline_version,
        author = 'Alex Griffing',
        author_email = 'argriffi@ncsu.edu',
        maintainer = 'Alex Griffing',
        maintainer_email = 'argriffi@ncsu.edu',
        url = 'http://github.com/argriffing/drosophiline',
        description = 'A collection of drosophila analysis scripts.',
        long_description = long_description,
        classifiers = classifiers,
        platforms = ['Linux'],
        packages = ['drosophiline'],
        scripts = scripts)
