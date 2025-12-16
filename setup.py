from setuptools import setup, find_packages

setup(
    name="kegg_crosstalk_analysis",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "networkx",
        "requests"
    ],
    description="KEGG pathway crosstalk analysis with Jaccard normalization",
    author="Durvi",
)
