import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="maliampi",
    version="0.5",
    author="Jonathan Golob",
    author_email="j-dev@golob.org",
    description="Maximum Likelihood Amplicon Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jgolob/maliampi",
    packages=setuptools.find_packages(),
    #dependency_links=[
    #    'https://github.com/jgolob/sciluigi/tarball/containertask#egg=sciluigi'
    #],
    install_requires=[
        'sciluigi @ git+git://github.com/jgolob/sciluigi@containertask#egg=sciluigui',
        'biopython',
        'pytz',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points= {
        'console_scripts': ['maliampi=maliampi.maliampi:main']
    }
)