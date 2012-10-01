from distutils.core import setup

# Find out what version is in this installation directory.
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), 'src/'))
import qecc as q
f=open('README.rst','r')

setup(
    name='QuaEC',
    version='{0}.{1}.{2}'.format(*q.__version__),
    url='http://cgranade.github.com/python-quaec/',
    author='Chris Granade and Ben Criger',
    author_email='cgranade@cgranade.com',
    package_dir={'': 'src'},
    packages=['qecc'],
    long_description=f.read()
)

f.close()