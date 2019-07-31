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

* Probability right to the end

    Using 16S rRNA gene amplicons to establish the structure and composition of a microbial community is inherently an uncertain process. This pipeline uses software and approaches that respect this uncertainty, and bring it right to the terminal analyses. 

* Runs the same locally (via docker), on slurm or PBS (via singularity), or on AWS.

    Maliampi's use of containerized software allows the pipeline to run the same on a variety of different computational resources. You have slurm, but your collaborators have access to a PBS cluster? No issue. On the go? Fine to run the pipeline's computationally heavier tasks via AWS.

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

### 0. Prepare our output and working directories
`mkdir -p working/ && mkdir -p output/sv/`

### 1. Create sequence variants
>Maliampi uses a manifest file to keep track of things. At a minimum this needs to have two columns: one for _specimen_ (a unique identifier) and _read__1_ (to point to the file path of the fastq file containing the forward read). Almost always there will be another column for the reverse read, _read__2_. 
> Another important column to include is _batch_, corresponding to a given extraction-and-sequencing run. 
> Take a look at manifest.csv to see this simple example of samples from the ATCC mock community.


```
maliampi --workers 2 \
sv_dada2 \
--manifest manifest.csv \
--working-dir working \
--destination-dir output/sv/
```

> A lot of log lines will spill by. Keep an eye out for 'xxx tasks pending' to get a sense of how much longer. Maliampi will run dada2 on each specimen in parallel, to the maximum number of workers you specified above at a time. 

If all succeeded, you should see
```
===== Luigi Execution Summary =====

Scheduled 70 tasks of which:
* 9 complete ones were encountered:
    - 9 LoadSpecimenReads(...)
* 61 ran successfully:
    - 9 BCCSpecimenReads(...)
    - 1 DADA2_Combine_Seqtabs(...)
    - 9 DADA2_DADA(...)
    - 9 DADA2_Dereplicate(...)
    - 9 DADA2_FilterAndTrim(...)
    ...

This progress looks :) because there were no failed tasks or missing dependencies

===== Luigi Execution Summary =====
```

> You can look to see the output files created, in 
`output\sv\`

> These include a fasta file, a matrix, and a map and weight file we will use later. 

### 2. Create a reference package
```
maliampi --workers 1 \
refpkg \
--working-dir working/ \
--refpkg-destdir output/ \
--refpkg-name tutorial \
--sequence-variants output/sv/dada2.sv.fasta \
--entrez-email your@email.com \
--repo-seq-info test.repo.seq_info.csv \
--repo-annotated-fasta test.repo.fasta
```

If this all works, you should see something like 
```
===== Luigi Execution Summary =====

Scheduled 11 tasks of which:
* 4 complete ones were encountered:
    - 1 AlignmentStoToFasta(...)
    - 1 CMAlignSeqs(...)
    - 1 CombineRepoMatches(...)
    - 1 ObtainCM(...)
* 7 ran successfully:
    - 1 BuildTaxtasticDB(...)
    - 1 CleanupTreeInfo(...)
    - 1 CombineRefpkg(...)
    - 1 ConfirmSeqInfoTaxonomy(...)
    - 1 RAxMLTree(...)
    ...

This progress looks :) because there were no failed tasks or missing dependencies

===== Luigi Execution Summary =====
```

The reference package will be packaged up at `output/refpkg/tutorial.date_time.refpkg.tgz` where date and time are stamps from when your reference package was made. 

### 3. Place SV on the Reference Package
Please note, you will need to replace date_time below with the timestamp from your reference package.

```
maliampi --workers 2 \
placement \
--working-dir working \
--destination-dir output \
--sequence-variants output/sv/dada2.sv.fasta \
--seq-map-csv output/sv/dada2.sv.map.csv \
--sv-weights-csv output/sv/dada2.sv.weights.csv \
--refpkg-tgz output/refpkg/tutorial.date_time.refpkg.tgz 
```

If all works, you should see something like
```
===== Luigi Execution Summary =====

Scheduled 16 tasks of which:
* 4 complete ones were encountered:
    - 1 LoadFastaSeqs(...)
    - 2 LoadFile(...)
    - 1 LoadRefpkgTGZ(...)
* 12 ran successfully:
    - 1 CMAlignSeqs(...)
    - 1 CombineAlignmentsSTO(...)
    - 1 ExtractRefpkgAlignment(...)
    - 1 Jplace_ADCL(...)
    - 1 Jplace_Alpha_Diversity(...)
    ...

This progress looks :) because there were no failed tasks or missing dependencies

===== Luigi Execution Summary =====
```

You can check out the placement outputs, which include a json file with the placement details, PCA, as well as alpha diversity and pairwise distance data, all in `output\placement`:
```
output/placement/
├── adcl.gz
├── alpha_diversity.csv
├── dedup.jplace
├── edpl.gz
├── kr_distance.csv
├── pca
│   ├── epca.proj
│   ├── epca.trans
│   ├── epca.xml
│   ├── lpca.proj
│   ├── lpca.trans
│   └── lpca.xml
└── redup.jplace.gz
```

### 4. Taxonomic classification
Please note, you will need to replace date_time below with the timestamp from your reference package.

```
maliampi --workers 1 \
classify \
--working-dir working \
--destination-dir output \
--sequence-variants output/sv/dada2.sv.fasta \
--seq-map-csv output/sv/dada2.sv.map.csv \
--sv-weights-csv output/sv/dada2.sv.weights.csv \
--jplace output/placement/dedup.jplace \
--refpkg-tgz output/refpkg/tutorial.date_time.refpkg.tgz
```

If all worked, you should see
```
INFO: 
===== Luigi Execution Summary =====

Scheduled 16 tasks of which:
* 5 complete ones were encountered:
    - 1 CombineAlignmentsSTO(...)
    - 3 LoadFile(...)
    - 1 LoadRefpkgTGZ(...)
* 11 ran successfully:
    - 6 GenerateTables(...)
    - 1 PlacementDB_AddSI(...)
    - 1 PlacementDB_Classify_SV(...)
    - 1 PlacementDB_MCC(...)
    - 1 PlacementDB_Prep(...)
    ...

This progress looks :) because there were no failed tasks or missing dependencies

===== Luigi Execution Summary =====
```

Look in the `output/classification/` directory to see what we've created:
> `placement.db` is an sqlite3 database.
```
├── placement.db
└── tables
    ├── by_specimen.class.csv
    ├── by_specimen.family.csv
    ├── by_specimen.genus.csv
    ├── by_specimen.order.csv
    ├── by_specimen.phylum.csv
    ├── by_specimen.species.csv
    ├── by_taxon.class.csv
    ├── by_taxon.family.csv
    ├── by_taxon.genus.csv
    ├── by_taxon.order.csv
    ├── by_taxon.phylum.csv
    ├── by_taxon.species.csv
    ├── tallies_wide.class.csv
    ├── tallies_wide.family.csv
    ├── tallies_wide.genus.csv
    ├── tallies_wide.order.csv
    ├── tallies_wide.phylum.csv
    └── tallies_wide.species.csv
```

