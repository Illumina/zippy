"""
Setup for Zippy (zippy)
"""
from setuptools import setup, find_packages

package = __import__('zippy')

def readme():
    with open('README.md') as f:
        return f.read()

try:
    long_desc = readme()
except IOError as err:
    long_desc = str(err)

try:
    version_str = open('version.txt').read()
except IOError as err:
    version_str = '2.0.0'

setup(
    name='zippy-pipeline',
    long_description=long_desc,
    url='https://github.com/Illumina/zippy',
    download_url='https://github.com/Illumina/zippy/archive/2.1.0.tar.gz',
    author='Aaron Wise',
    author_email='quejebo@gmail.com',
    license='Apache 2.0',
    install_requires=[
        'distribute',
        'commentjson',
        'ast',
        'pybloomfiltermmap',
        'pyflow'
        ],
    dependency_links=['https://github.com/Illumina/pyflow/releases/download/v1.1.17/pyflow-1.1.17.tar.gz'],
    version=version_str,
    description='NGS Pipeline Prototyping Tool for Python',
    include_package_data=True,
    packages=find_packages(),
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=[],
    requires=[],
    zip_safe=False)
