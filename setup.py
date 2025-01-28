from setuptools import setup, find_packages
from src import __version__ as version

setup(
    name='mhcbooster',
    version=str(version),
    # url='https://github.com/caronlab/mhcbooster',
    license='',
    author='Ruimin Wang',
    author_email='',
    description='',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    entry_points={
        'console_scripts': ['mhcbooster = src.interface.mhcbooster_cli:run',
                            'mhcbooster-gui = src.interface.mhcbooster_gui']
    },
    packages=find_packages(),
    python_requires='==3.10',
    install_requires=['mhcflurry>=2.0.2', 'mhcnames', 'tensorflow==2.16.1', 'scikit-learn', 'pandas', 'numpy',
                      'pandas', 'tqdm', 'pyteomics', 'matplotlib', 'lxml', 'tensorflow-probability',
                      'xmltodict', 'im2deep', 'peptdeep']
)
