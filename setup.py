from setuptools import setup, find_packages, find_namespace_packages
from src import __version__ as version

def read_requirements(filename="requirements.txt"):
    with open(filename, "r") as f:
        return [line.strip() for line in f.readlines() if line.strip() and not line.startswith("#")]

setup(
    name='mhcbooster',
    version=str(version),
    description='',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url='https://github.com/caronlab/mhcbooster',
    author='Ruimin Wang',
    author_email='ruimin.wang@yale.edu',
    entry_points={
        'console_scripts': ['mhcbooster = src.interface.mhcbooster_cli:run',
                            'mhcbooster-gui = src.interface.mhcbooster_gui:run',
                            'mhcbooster-package-installer = src.utils.package_installer:install']
    },
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    packages=['src', 'src.adapter', 'src.interface', 'src.model', 'src.predictors', 'src.report', 'src.utils'],
    python_requires='==3.10',
    install_requires=read_requirements(),
    include_package_data=True,
    license='GPL-3.0'
)
