# MaLiAmPi
### Maximum Likelihood Amplicon Pipeline: An amplicon (PCR / 16S) microbiome pipeline.

## Introduction
Maliampi is a phylogenetic placement-based pipeline for handling 16S rRNA amplicons (the most common type of microbiome study data).

Maliampi's overall philosophy is to embrace and work with the inherent ambiguity involved in any PCR / seqencing based approach, as well as rationally deal with the limits of the available references mapping between sequence-space and taxonomy. 

### A typical workflow:
1) Create sequence variants
2) Make a reference package
3) Place on the reference package
4) Classify

### Advantages of Maliampi
* Containerized dependencies

    Maliampi uses docker containers for all of the software called by the pipeline. This means no fussing with getting archane academic software matched up properly with dependencies, or issues with reproducing the pipeline environment on your system.

* Modular design
    
    While maliampi can start with fastq files straight from the sequencer, you can also opt to use your own sequence variant / OTU generation strategy and just use maliampi to make a reference package. 



## Installation
Right now, maliampi relies upon a [custom-fork](https://github.com/jgolob/sciluigi/tree/containertask) of [sciluigi](https://github.com/pharmbio/sciluigi). That means one cannot simply install with pip. (Pull request is pending).

### Dependencies
* Python 3
* Docker (if you want to run locally)

## Installation
> For a variety of reasons, maliampi works best on an unix-like operating system (linux, MacOS, [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)). Directions here are based on an unix-like system. 

1) Install python3 with pip.

    Ubuntu: `sudo apt-get install python3-pip`

    Mac: I suggest via [brew](https://docs.python-guide.org/starting/install3/osx/).

2) (Optional) Install [docker engine community edition](https://hub.docker.com/search/?type=edition&offering=community).

    This is needed if you want to run things locally for testing and is strongly suggested.

    Windows + Docker is not a happy marriage right now--particularly when attempting to bind local directories in containers. Thus, maliampi does not work well on windows, even with docker and python installed. 

3) Create a python virtual environment (strongly recommended).

    i) `python3 -m venv maliampi-env`

    ii) `source maliampi-env/bin/activate`

    iii) `pip install pip --upgrade`


4) Obtain and install the custom fork of sciluigi

    i. [download sciluigi as a zip](https://github.com/jgolob/sciluigi/archive/containertask.zip)

    ii. Unzip: 
    
    `unzip containertask.zip`
    
    iii. Install using pip: 
    
    `pip3 install sciluigi-containertask/`

    iv. Copy over and edit the containerinfo.ini file.
    
    `mkdir ~/.sciluigi && cp sciluigi-containertask/example-config/containerinfo.ini ~/.sciluigi/`

    v. Cleanup `rm -r sciluigi-containertask && rm containertask.zip`

5) Obtain and install maliampi

    i. [download maliampi as a zip](https://github.com/jgolob/maliampi/archive/master.zip)

    ii. Unzip: 
    
    `unzip master.zip`
    
    iii. Install using pip: 
    
    `pip3 install maliampi-master/`

    iv. Move the practice data

    `mv maliampi-master/maliampi-practice-data/ .`

    v. Cleanup 
    
    `rm -r maliampi-master && rm master.zip`

## Tutorial
>The tutorial presumes you have maliampi installed as above (and if in a virtual environment, the virtual environment activated).

> The tutorial uses the data in the maliampi-practice-data folder from the repository, with the command lines below assuming you are within that directory.

### I. Create sequence variants
>Maliampi uses a manifest file to keep track of things. At a minimum this needs to have two columns: one for _specimen_ (a unique identifier) and _read__1_ (to point to the file path of the fastq file containing the forward read). Almost always there will be another column for the reverse read, _read__2_. Take a look at manifest.csv to see this simple example of samples from the ATCC mock community.

1. Create our working and output directories

`mkdir -p working/ && mkdir -p output/sv/

2. Run maliampi in sv_dada2 mode:
    ```
    maliampi --workers 2 \
    sv_dada2 \
    --manifest manifest.csv \
    --working-dir working \
    --destination-dir output/sv/
    ```

### II. Create a reference package
    ```
    maliampi --workers 1 \
    refpkg \
    --manifest manifest.csv \
    --working-dir working \
    --destination-dir output/sv/
    ```
