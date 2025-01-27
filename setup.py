from setuptools import setup, find_packages
from mhcvalidator import __version__ as version

setup(
    name='mhcbooster',
    version=str(version),
    packages=find_packages(),
    # url='https://github.com/caronlab/mhc-validator',
    license='',
    author='Ruimin Wang',
    author_email='',
    description='',
    install_requires=['mhcflurry>=2.0.2', 'mhcnames', 'tensorflow==2.16.1', 'scikit-learn', 'pandas', 'numpy',
                      'pandas', 'tqdm', 'pyteomics', 'matplotlib', 'lxml', 'tensorflow-probability',
                      'xmltodict', 'im2deep', 'peptdeep'],
    entry_points={
        'console_scripts': ['mhcvalidator = mhcvalidator.command_line:run']
    },
    python_requires='==3.10'
)
