from distutils.command.sdist import sdist as sdist_orig
from distutils.errors import DistutilsExecError
from setuptools import setup, find_packages
from os import path

PKG_NAME = "Cas13gRNAtor"
MOD_NAME = "cas13gRNAtor"

DESCRIPTION = """ 
Cas13gRNAtor
RfxCas13d guide RNA Scoring and Selection
Based on Guide Efficacy, Shannon Entropy and Conservation Scores of gRNA 
"""

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as fh, open(path.join(here, "requirements.txt")) as req:
    install_requires = [pkg.strip() for pkg in req]

__version__ = ""
exec(open("{}/_version.py".format(MOD_NAME)).read())

data = [
    "misc/bowtie",
    "misc/bowtie-build",
    "misc/RfxCas13d_GuideScoring.R",
    "misc/RNAfold",
    "misc/RNAhybrid",
    "misc/RNAplfold",
    "misc/RNAhyb.sh",
    "misc/data/Cas13designGuidePredictorInput.csv",
    "misc/data/LocalNTdensityCorrelations.txt",
    "misc/data/test.fa",
]

#FASTA = f"{MOD_NAME}/misc/data/MN908947_NY1-PV08001.S.fasta"
#os.system(f"Rscript {MOD_NAME}/{misc}/RfxCas13d_GuideScoring.R {fasta} {MOD_NAME}/{misc}/data/Cas13designGuidePredictorInput.csv true")
'''
class sdist(sdist_orig):
    def run(self):
        try:
            self.spawn(['Rscript', f'{MOD_NAME}/misc/RfxCas13d_GuideScoring.R', FASTA, f"{MOD_NAME}/misc/data/Cas13designGuidePredictorInput.csv", "true"])
        except DistutilsExecError:
            self.warn('Command to install R packages failed')
        super().run()

    #cmdclass={'install': sdist},
'''
setup(
    name=PKG_NAME,
    version=__version__,
    author="Muhammad Irfan",
    author_email="muhammad_irfan@gis.astar.edu.sg",
    description="RfxCas13d guide RNA Scoring and Selection",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/chewlabSB2",
    packages=find_packages(exclude=["*.test", "*.test.*", "test.*", "test"]),
    package_dir={'cas13gRNAtor': 'cas13gRNAtor'}, 
    package_data={'cas13gRNAtor': data},
    entry_points={
        "console_scripts": [
            "cas13gRNAtor={}.cas13gRNAtor:main".format(MOD_NAME),
            "cas13gRNAscore={}.cas13gRNAscore:main".format(MOD_NAME),
            "cas13gRNAvalidate={}.cas13gRNAvalidate:main".format(MOD_NAME),
            "cas13gRNAinstall={}.installRpackage:main".format(MOD_NAME),
            "cas13gRNApanvirus={}.cas13gRNApanvirus:main".format(MOD_NAME),
            #"cas13gRNApredict={}.gRNA_score:main".format(MOD_NAME),
        ],
    },
    
    install_requires=install_requires,
    include_package_data=True,
    python_requires=">=3.5",
)

import os 
os.system("cas13gRNAinstall")