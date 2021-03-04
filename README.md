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

* Probability right to the end

    Using 16S rRNA gene amplicons to establish the structure and composition of a microbial community is inherently an uncertain process. This pipeline uses software and approaches that respect this uncertainty, and bring it right to the terminal analyses.

    Alpha-diversity, phylogenetic distance, and taxonomic classification are all done respecting the uncertainty involved in each step.

* Runs the same locally (via docker), on slurm or PBS (via singularity), or on the cloud.

    Maliampi's use of containerized software and nextflow allows the pipeline to run the same on a variety of different computational resources. You have slurm, but your collaborators have access to a PBS cluster? No issue. On the go? Fine to run the pipeline's computationally heavier tasks via AWS while running simpler steps on docker right on your computer.

* Modular design
    
    (TODO)
    While maliampi can start with fastq files straight from the sequencer, you can also opt to use your own sequence variant / OTU generation strategy and just use maliampi to make a reference package. 

## Installation
Maliampi uses [nextflow](https://www.nextflow.io/) as a workflow engine. 

Using maliampi requires nextflow to be installed and configured on your computer. Here we will go through the steps of doing so.

(A previous attempt at maliampi used sciluigi. This is now deprecated, with the code remaining around for now.)

### Dependencies
* [nextflow](https://www.nextflow.io/), which in turn requires the java runtime environment.
* Docker (if you want to run locally)

## Installation
> For a variety of reasons, maliampi works best on an unix-like operating system (linux, MacOS, [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)). Directions here are based on an unix-like system. 

1) Install nextflow

    [Nextflow installation directions](https://www.nextflow.io/index.html#GetStarted) can be followed and are only three steps:
    
    i) Install the java runtime environment (at least version 8). I suggest [OpenJDK 8](https://adoptopenjdk.net/).

    ii) Use curl to download the latest nextflow binary:
    
    `curl -s https://get.nextflow.io | bash`

    iii) Test nextflow

    `./nextflow run hello`

2) (Optional) Install [docker engine community edition](https://hub.docker.com/search/?type=edition&offering=community).

    This is needed if you want to run things locally for testing and is strongly suggested. Directions can be found on the docker site.

    >Windows + Docker is not a happy marriage right now--particularly when attempting to bind local directories in containers. 

3) Create (or edit) a common nextflow config file

    We need to tell nextflow what sort of computational resources are available in your setup.

    Different steps in maliampi require different resources. Maliampi gives each step (or process) a label:

    - 'io_limited': processes that are not multithreaded and are generally limited by the speed files can be read and written.

    - 'io_limited_local': processes that are I/O limited and should also be executed locally.

    - 'io_mem': processes that are limited by the speed of reading and writing files, but also require a lot of memory.

    - 'mutlithread': Tasks that have relatively modest cpu and memory requirements, but can benefit from multiple cpus and are not limited by reading and writing files.

    - 'mem_veryhigh': Multithreaded tasks that require very large amounts of memory.

    Maliampi works best when the number of cpus and memory for each of these labels is established in the common nextflow config file. 
    
    In addition the config file should tell nextflow to use docker locally (or AWS / GCS / slurm / PBS / singularity depending on your setup). Details of this can be found in the nextflow documentation. You can seek out the help of your local high-performance computing staff and the nextflow documentation.

    An example config file that uses docker locally (and can be copied into `~/.nextflow/config` to get you started):
```groovy
profiles{
    standard {
        process {
            executor = 'local'
            withLabel: 'io_limited' {
                cpus = 1
                memory = 4.GB
            }
            withLabel: 'io_limited_local' {
                cpus = 1
                memory = 4.GB
            }
            withLabel: 'io_mem' {
                cpus = 1
                memory = 16.GB
            }
            withLabel: 'multithread' {
                cpus = 4
                memory = 16.GB
            }
            withLabel: 'mem_veryhigh' {
                cpus = 2
                memory = 32.GB
            }
        }
        docker {
            enabled = true
        }
    } // end standard profile
    hybrid {  // start hybrid profile combining local and awsbatch
        process {
            withLabel: 'io_limited' {
                cpus = 1
                memory = 4.GB
                executor = 'awsbatch'
            }
            withLabel: 'io_limited_local' {
                cpus = 1
                memory = 4.GB
                executor = 'local'
                docker.enabled = true
            }
            withLabel: 'io_mem' {
                cpus = 1
                memory = 16.GB
                executor = 'awsbatch'
            }
            withLabel: 'multithread' {
                cpus = 4
                memory = 16.GB
                executor = 'awsbatch'
            }
            withLabel: 'mem_veryhigh' {
                cpus = 2
                memory = 32.GB
                executor = 'awsbatch'
            }
        }
    } // end hybrid profile
}
```

   >If your computer setup is not that rich in memory, feel free to cut down any of the specific memory labels.
   
   >This setup is useful for testing and work with modest data sets, including our test data. For larger data, setting up nextflow to use a high-performance computing cluster and/or cloud computing (AWS, GCS, etc) will be required in most situations.


## Tutorial
>The tutorial presumes you have nextflow installed as above along with the `~/.nextflow/config` file set up as described above.


### 1. Obtain the practice data
>The maliampi repository contains some sample data we can work with. We will download the current version of the repository to obtain it, save the practice data, then cleanup the rest.
```
wget https://github.com/jgolob/maliampi/archive/master.zip
unzip master.zip
mv maliampi-master/maliampi-practice-data/ .
rm -r maliampi-master && rm master.zip
cd maliampi-practice-data
```

### 2. Prepare our output and working directories
```
mkdir -p working/ && mkdir -p output/sv/
```

### 3. (Optional) Look at the manifest file
>Maliampi uses a manifest file to keep track of things. At a minimum this needs to have two columns: one for _specimen_ (a unique identifier) and _R1_ (to point to the file path of the fastq file containing the forward read). Almost always there will be another column for the reverse read, _R2_. 

> Another important column to include is _batch_, corresponding to a given extraction-and-sequencing run. 

> Take a look at manifest.csv to see this simple example of samples from the ATCC mock community.


### 4. Use `nextflow` to show the maliampi command-line options

```
nextflow run jgolob/maliampi -latest --help
```

You should see:
```
N E X T F L O W  ~  version 19.07.0
Pulling jgolob/maliampi ...
downloaded from https://github.com/jgolob/maliampi.git
Launching `jgolob/maliampi` [modest_saha] - revision: fc96dbdbfb [master]
Usage:

nextflow run jgolob/maliampi <ARGUMENTS>

Required Arguments:
    --manifest            CSV file listing samples
                            At a minimum must have columns:
                                specimen: A unique identifier 
                                read__1: forward read
                                read__2: reverse read fq

                            optional columns:
                                batch: sequencing / library batch. Should be filename safe
                                index__1: forward index file (for checking demultiplexing)
                                index__2: reverse index file
    --repo_fasta          Repository of 16S rRNA genes.
    --repo_si             Information about the 16S rRNA genes.
    --email               Email (for NCBI)
Options:
  Common to all:
    --output              Directory to place outputs (default invocation dir)
                            Maliampi will create a directory structure under this directory
    -w                    Working directory. Defaults to `./work`
    -resume                 Attempt to restart from a prior run, only completely changed steps

  SV-DADA2 options:
    --trimLeft              How far to trim on the left (default = 0)
    --maxN                  (default = 0)
    --maxEE                 (default = Inf)
    --truncLenF             (default = 0)
    --truncLenR             (default = 0)
    --truncQ                (default = 2)

  Ref Package required:
    --repo_fasta            FASTA file containing reference sequences (required)
    --repo_si               CSV file with information about the repo reads (required)
  Ref Package options (defaults generally fine):
    --repo_min_id           Minimum percent ID to a SV to be recruited (default = 0.8)
    --repo_max_accepts      Maximum number of recruits per SV (default = 10)
    --cmalign_mxsize        Infernal cmalign mxsize (default = 8196)
    --raxml_model           RAxML model for tree formation (default = 'GTRGAMMA')
    --raxml_parsiomony_seed (default = 12345)
    --taxdmp                Path to taxdmp.zip. If not provided, it will be downloaded

  Placement / Classification Options (defaults generally fine):
    --pp_classifer                  pplacer classifer (default = 'hybrid2')
    --pp_likelihood_cutoff          (default = 0.9)
    --pp_bayes_cutoff               (default = 1.0)
    --pp_multiclass_min             (default = 0.2)
    --pp_bootstrap_cutoff           (default = 0.8)
    --pp_bootstrap_extension_cutoff (default = 0.4)
    --pp_nbc_boot                   (default = 100)
    --pp_nbc_target_rank            (default = 'genus')
    --pp_nbc_word_length            (default = 8)
    --pp_seed                       (default = 1)
```

> maliampi takes a lot of options, but only a few are required:

> --manifest: the csv file telling us where to find files

> --repo_fasta: A fasta file of 16S rRNA genes

> --repo_si: A CSV file linking the 16S rRNA gene IDs in the fasta file with the NCBI taxonomy.

> --email: a valid email address to use with NCBI

> nextflow will automatically download maliampi from the repository on github. `-latest` ensures the latest version is always downloaded.


### 5. Use `nextflow` to run maliampi on the practice data
```
nextflow run jgolob/maliampi -latest -resume \
--manifest manifest.csv \
--output tutorial_output \
--email your@email.com \
--repo_fasta test.repo.fasta \
--repo_si test.repo.seq_info.csv 
```
Maliampi will run. It will take a bit to complete even on the test data.

If this all works, you should see something like:
```
N E X T F L O W  ~  version 19.07.0
Pulling jgolob/maliampi ...
Already-up-to-date
Launching `jgolob/maliampi` [furious_lalande] - revision: 2160a88d19 [master]
executor >  local (96)
[03/d23f0b] process > barcodecop (9)                 [100%] 9 of 9 ✔
[e0/e314a2] process > dada2_ft (9)                   [100%] 9 of 9 ✔
[1e/9fc05c] process > dada2_derep (9)                [100%] 9 of 9 ✔
[9c/2cbe37] process > dada2_learn_error (3)          [100%] 6 of 6 ✔
[1b/e8a4b9] process > dada2_derep_batches (6)        [100%] 6 of 6 ✔
[de/8cbc57] process > dada2_dada (6)                 [100%] 6 of 6 ✔
[2e/5034c2] process > dada2_demultiplex_dada (6)     [100%] 6 of 6 ✔
[62/fa4431] process > dada2_merge (9)                [100%] 9 of 9 ✔
[7e/2748f5] process > dada2_seqtab_sp (9)            [100%] 9 of 9 ✔
[67/003c58] process > dada2_seqtab_batch_combine (2) [100%] 3 of 3 ✔
[de/6a7a43] process > dada2_seqtab_combine_all       [100%] 1 of 1 ✔
[e2/90044f] process > dada2_remove_bimera            [100%] 1 of 1 ✔
[0a/920e9f] process > dada2_convert_output           [100%] 1 of 1 ✔
[bd/3bd963] process > output_failed                  [100%] 1 of 1 ✔
[7b/68d107] process > refpkgSearchRepo               [100%] 1 of 1 ✔
[cb/c47830] process > filterSeqInfo                  [100%] 1 of 1 ✔
[7b/5fdc72] process > DlBuildTaxtasticDB             [100%] 1 of 1 ✔
[6f/da64ae] process > confirmSI                      [100%] 1 of 1 ✔
[ba/01d471] process > alignRepoRecruits              [100%] 1 of 1 ✔
[de/b68cf4] process > convertAlnToFasta              [100%] 1 of 1 ✔
[f7/7251c8] process > raxmlTree                      [100%] 1 of 1 ✔
[d3/9099c7] process > raxmlTree_cleanupInfo          [100%] 1 of 1 ✔
[d0/657abb] process > taxtableForSI                  [100%] 1 of 1 ✔
[2b/9ddac1] process > obtainCM                       [100%] 1 of 1 ✔
[d6/eabb78] process > combineRefpkg                  [100%] 1 of 1 ✔
[0e/075e7e] process > alignSV                        [100%] 1 of 1 ✔
[5d/be1074] process > extractRefpkgAln               [100%] 1 of 1 ✔
[ed/97c83b] process > combineAln_SV_refpkg           [100%] 1 of 1 ✔
[11/4b6835] process > pplacerPlacement               [100%] 1 of 1 ✔
[57/5bf9d9] process > pplacerReduplicate             [100%] 1 of 1 ✔
[91/afe8bd] process > pplacerADCL                    [100%] 1 of 1 ✔
[48/8e6d1b] process > pplacerEDPL                    [100%] 1 of 1 ✔
[25/5d3693] process > pplacerPCA                     [100%] 1 of 1 ✔
[90/ee0880] process > pplacerAlphaDiversity          [100%] 1 of 1 ✔
[3a/348adb] process > pplacerKR                      [100%] 1 of 1 ✔
[34/e44686] process > classifyDB_Prep                [100%] 1 of 1 ✔
[14/afd81c] process > classifySV                     [100%] 1 of 1 ✔
[d8/0d2f39] process > classifyMCC                    [100%] 1 of 1 ✔
[d7/d5ca1e] process > classifyTables (6)             [100%] 6 of 6 ✔
Completed at: 12-Sep-2019 14:58:56
Duration    : 13m 16s
CPU hours   : 0.9
Succeeded   : 96
Cached      : 0
```

### 6. Look at the output
There are four subdirectories in the output:
```
tutorial_output/
├── classify
├── placement
├── refpkg
└── sv
```

### 6A. `sv` for sequence variants
```
tutorial_output/sv/
├── dada2.combined.seqtabs.nochimera.csv
├── dada2.combined.seqtabs.nochimera.rds
├── dada2.specimen.sv.long.csv
├── dada2.sv.fasta
├── dada2.sv.map.csv
├── dada2.sv.shared.txt
├── dada2.sv.weights.csv
├── errM
│   ├── m22_p4
│   │   ├── m22_p4.R1.errM.csv
│   │   ├── m22_p4.R1.errM.rds
│   │   ├── m22_p4.R2.errM.csv
│   │   └── m22_p4.R2.errM.rds
│   ├── m23_p3
│   │   ├── m23_p3.R1.errM.csv
│   │   ├── m23_p3.R1.errM.rds
│   │   ├── m23_p3.R2.errM.csv
│   │   └── m23_p3.R2.errM.rds
│   └── m48_p4
│       ├── m48_p4.R1.errM.csv
│       ├── m48_p4.R1.errM.rds
│       ├── m48_p4.R2.errM.csv
│       └── m48_p4.R2.errM.rds
└── failed_specimens.csv
```
Some highlights here:
- `dada2.sv.fasta` contains the sequence variants in FASTA format.
- `dada2.specimen.sv.long.csv` contains the SV count by community in long format (i.e. a CSV file with columns community, sv, and count).
- `dada2.sv.shared.txt` contains the SV counts in a TSV format similar to `mothur` sharefile. 
- `failed_specimens.csv` contains all of the specimens that failed at some point (files empty, none survive filter-trim, merging, etc). The reason is listed as well. This is useful for troubleshooting.
- The `errM` subdirectory contains the error models for each batch. 
> By default, maliampi will try to ignore failed specimens and finish with as many specimens as possible.

### 6B. `refpkg` for the reference package
```
tutorial_output/refpkg/
└── refpkg.tgz
```
A reference package is used for placement, and stored in tar-gzip format. Within the tar-gzip file is:
```
drwxr-xr-x  0 root   root        0 Sep 12 14:57 ./
-rw-r--r--  0 root   root     2217 Sep 12 14:57 ./CONTENTS.json
-rw-r--r--  0 root   root      381 Sep 12 14:57 ./phylo_modeldqujt6x8.json
-rw-r--r--  0 root   root     7402 Sep 12 14:57 ./RAxML_bestTree.refpkg
-rw-r--r--  0 root   root     2239 Sep 12 14:57 ./RAxML_info.refpkg
-rw-r--r--  0 root   root   207936 Sep 12 14:57 ./recruits.aln.fasta
-rw-r--r--  0 root   root   252268 Sep 12 14:57 ./recruits.aln.sto
-rw-r--r--  0 root   root  1012905 Sep 12 14:57 ./refpkg.cm
-rw-r--r--  0 root   root    28026 Sep 12 14:57 ./refpkg.seq_info.corr.csv
-rw-r--r--  0 root   root    27359 Sep 12 14:57 ./refpkg.taxtable.csv
-rw-r--r--  0 root   root     4340 Sep 12 14:57 ./treejkkusc5k.tre
```
- `CONTENTS.json` a manifest.
- `RAxML_bestTree.refpkg` a tree in newick format.
- `recruits.aln.fasta` an alignment of the full-length 16s rRNA gene recruits (FASTA format).
- `recruits.aln.sto` an alignment of the full-length 16s rRNA gene recruits (Stockholm format).
- `refpkg.cm` the covariance matrix used for the alignment.
- `refpkg.seq_info.corr.csv` A sequence information file telling us about the full length rRNA genes (inlcuding their taxonomic assignment)
- `refpkg.taxtable.csv` a mothur-style table of taxonomic information for each full-length 16s rRNA gene in the reference.

### 6C. `placement` for the phylogenetic placement of sequence variants

```
tutorial_output/placement/
├── adcl.csv.gz
├── alpha_diversity.csv.gz
├── dedup.jplace
├── edpl.csv.gz
├── kr_distance.csv.gz
├── pca
│   ├── epca.proj
│   ├── epca.trans
│   ├── epca.xml
│   ├── lpca.proj
│   ├── lpca.trans
│   └── lpca.xml
└── redup.jplace.gz
```
- `dedup.jplace` is the deduplicated placement file in JSON format.
- `redup.jplace.gz` is the reduplicated placement file.
- `alpha_diversity.csv.gz` has alpha diversity information for each community.
- `kr_distance.csv.gz` has the weighted pairwise phylogenetic distance between communities in long format.
- The `pca\` subdirectory has the PCA coordinates for each community.
- `adcl.csv.gz` and `edpl.csv.gz` contains information about each sequence variant's placement (Average distance to the closest leaf and pendant length respectively). High ADCL for a SV is indicative of poor representation in the reference package.

### 6D. `classify` for the taxonomic classifications of sequence variants and for each specimen / community
```
tutorial_output/classify/
├── classify.mcc.db
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
- `classify.mcc.db` contains all of the classification information in a sqlite3 database. The most important tables are `multiclass` and `multiclass_concat`. The documentation for `pplacer` contains further information.
- the `tables/` subdirectory contains by-community and by target-rank classification information. 
- `tables/tallies_wide.*` contains the by-community classified output in wide format.
- `tables/by_specimen.*` contains the by-community taxon counts and fractional abundance in long format. 