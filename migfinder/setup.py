from setuptools import setup

setup(name='migfinder',
      version='1.0',
      description='Metagenomic Integron-associated Gene finder',
      url='http://github.com/maribuon/MIG-finder',
      author='Mariana Buongermino Pereira',
      author_email='maribuon@gmail.com',
      license='BSD-3-Clause',
      install_requires=['biopython==1.69'],
      zip_safe=False,
			package=['migfinder'],
			package_dir={'migfinder': 'migfinder/'},
			package_data={"migfinder": ["cm_model/selection109_oriR.cm"]
			},
)
