from setuptools import setup, find_packages

setup(
    name="timoris",
    version="0.1.0",
    author="Shu Zhou",
    author_email="shuzh@umich.edu",
    description="A tool for detecting tubule-like structures in spatial transcriptomics data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/shuzh/TIMORIS",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "matplotlib",
        "seaborn",
        "scikit-learn",
        "scanpy",
        "napari",
        "opencv-python",
        "pyvista"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'timoris-detect=timoris.detect_tubules:main',
            'timoris-visualize=timoris.visualize_results:main'
        ],
    },
    include_package_data=True,
)
