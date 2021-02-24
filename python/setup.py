import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "vcfLoader",
    version = "0.1",
    author = "Carles Hernandez-Ferrer and Davide Piscia",
    author_email = "davide.piscia@cnag.crg.eu",
    description = (""),
    license = "MIT",
    keywords = ["API", "VCF", "gVCF", "genome", "annotation"],
    url = "",
    packages=['rdconnect'],
    long_description=read('README')
)
