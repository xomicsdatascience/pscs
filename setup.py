import setuptools


def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()

setuptools.setup(
    name="pscs",
    version="0.0.4",
    author="Meyer Lab",
    author_email="",
    description="",
    long_description="",
    url="",
    classifiers=[
        "Development Status :: 2 - Alpha",
        'Intended Audience :: Science/Research',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires='>=3.10',
    license="MIT",
    install_requires=[
        'pyteomics>=4.4.1',
        'matplotlib<3.7',
        'numba>=0.53.1',
        'numpy>=1.20.1',
        'pandas>=1.2.2',
        'Bio>=0.4.1',
        'PyQt5>=5.15.4',
        'lxml>=4.6.2',
        'flask',
        'anndata',
        'scanpy',
        'plotly'
    ],
)
