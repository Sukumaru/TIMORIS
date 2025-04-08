from setuptools import setup, find_packages

setup(
    name="timoris",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scanpy",
        "matplotlib",
        "scikit-learn",
        "networkx",
        "scipy"
    ],
    author="Your Name",
    description="Tubule inference from spatial transcriptomics using ray-based centroid detection and refinement",
    license="MIT",
)