from distutils.core import setup, Command
import sys
from geo2fastq.config import VERSION

DESCRIPTION = """
Download GEO sequencing experiments and process to map and create track hubs.
"""

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass
        
    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)


setup(name = 'geo2fastq',
	version = VERSION,
	description = DESCRIPTION,
	author='Simon van Heeringen',
	author_email='s.vanheeringen@ncmls.ru.nl',
	license='MIT',
	packages=[
		'geo2fastq'
	],
	scripts=[
		"scripts/geo2trackhub",
	],
     requires=[], #TODO list required packages here
     platforms=["linux"],
     data_files=[
        ('config', ["config/geo2fastq.yaml"]),
     ],
     #include_package_data=True,
     tests_require=['pytest'],
     cmdclass = {'test': PyTest},
)
