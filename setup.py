from setuptools import setup

setup(
    name="fragment_analyzer",
    version="0.1.0",
    description="Fragment Analysis package in python!",
    url="https://github.com/Clinical-Genomics-Umea/ladder_map",
    author="William Rosenbaum and Pär Larsson",
    author_email="william.rosenbaum@umu.se",
    license="MIT",
    packages=["fragment_analyzer"],
    install_requires=[
        "pandas",
        "numpy",
        "scikit-learn",
        "matplotlib",
        "networkx",
        "lmfit",
        "scipy",
        "biopython",
    ],
)
