import glob
from setuptools import setup, find_packages

VERSION = '0.1'
DESCRIPTION = 'MIBiG conversion from 1.4 to latest schema'

def find_data_files():
    data_files = []
    for pathname in glob.glob("antismash/**/*", recursive=True):
        if pathname.endswith('.pyc'):
            continue
        if pathname.endswith('.py'):
            continue
        if '__pycache__' in pathname:
            continue
        data_files.append(pathname)
    return data_files

setup(
    name="mibig-converter",
    version="0.1",
    author="",
    author_email="",
    description=DESCRIPTION,
    entry_points={
        'console_scripts': [
            'convert-mibig=converter.__main__:entrypoint',
        ],
    },
    package_data={
        'converter': find_data_files(),
    },
    long_description="",
    packages=find_packages(),
    install_requires=[
        "jsonschema",
    ],
)
