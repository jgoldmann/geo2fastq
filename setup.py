from setuptools import setup
from geo2fastq.config import VERSION

DESCRIPTION = """
Download GEO sequencing experiments and process to map and create track hubs.
"""

setup(name = 'geo2fastq',
	version = VERSION,
     description = DESCRIPTION,
     long_description=open('README.md').read(),
	author='Simon van Heeringen',
	author_email='s.vanheeringen@ncmls.ru.nl',
      maintainer='Jakob M Goldmann',
      maintainer_email='jakob.goldmann@mail.de',
	license='MIT',
	packages=[
		'geo2fastq'
	],
	scripts=[
		"scripts/geo2trackhub",
           "scripts/geo2fastq",
	],
	install_requires=["biopython", 'PyYaml', 'pp'], 
	platforms=["linux"],
	data_files=[
	('config', ["config/geo2fastq.yaml"]),
	],
     test_suite='nose.collector',
     tests_require=['nose'],
	)
