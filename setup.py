from setuptools import setup

setup(name='migfinder',
      version='1.1',
      description='Metagenomic Integron-associated Gene finder',
      url='http://github.com/maribuon/MIG-finder',
      author='Mariana Buongermino Pereira',
      author_email='maribuon@gmail.com',
      license='BSD-3-Clause',
      install_requires=[
        'importlib_resources; python_version < "3.9"',
        'biopython>=1.79',
      ],
      zip_safe=False,
      package=['migfinder'],
      package_dir={'migfinder': 'migfinder/'},
      package_data={"migfinder": ["cm_model/selection109_oriR.cm"]
      },
      include_package_data=True,
      entry_points={
        'console_scripts': [
          'migfinder=migfinder.migfinder:migfinder_cli',
        ],
      },
)
