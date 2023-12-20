# boostDM-CH pipeline

**BoostDM-CH** is a method to score all possible single base substitutions in clonal hematopoiesis (CH) driver genes for their potential to drive CH.

The method's rationale and analysis results thereof are described in this study:

  *Identification of Clonal Hematopoiesis Driver Mutations through In Silico Saturation Mutagenesis*  
  Santiago Demajo, Joan Enric Ramis-Zaldivar, Ferran Muiños, Miguel L Grau, Maria Andrianova, Nuria Lopez-Bigas, Abel Gonzalez-Perez  
  URL: [https://doi.org/10.1101/2023.12.13.23299893](https://doi.org/10.1101/2023.12.13.23299893)

The method is based on two previous projects: **IntOGen-CH** and **boostDM**

**IntOGen-CH** pipeline provides a compendium of CH driver genes and mutational features upon using the **IntOGen** pipeline with a catalog of high VAF somatic mutations called in normal blood from tumor samples by reverse variant calling, which has been described in this study:

  *Discovering the drivers of clonal hematopoiesis*  
  Oriol Pich, Iker Reyes-Salazar, Abel Gonzalez-Perez, Nuria Lopez-Bigas  
  URL: [https://doi.org/10.1038/s41467-022-31878-0](https://doi.org/10.1038/s41467-022-31878-0)

**boostDM** is a methodology to annotate mutations in cancer genes for their potential to drive tumorigenesis, which has been described in another study:

  *In silico saturation mutagenesis of cancer genes*  
  Ferran Muiños, Francisco Martinez-Jimenez, Oriol Pich, Abel Gonzalez-Perez, Nuria Lopez-Bigas    
  URL: [https://www.nature.com/articles/s41586-021-03771-1](https://www.nature.com/articles/s41586-021-03771-1)

## Content
This repo contains the source code to reproduce the training and prediction steps of the boostDM-CH pipeline,
starting from the output data coming after the pipeline of IntOGen-CH

## Prerequisites

#### HPC

It is strongly recommended to run this pipeline in an HPC environment.

#### Singularity

The pipeline makes use of Singularity to handle software dependencies. 

* Create a singularity image using the [Singularity](https://github.com/bbglab/boostdm-pipeline/blob/ch/containers_build/boostdm/Singularity) recipe:
```
$ cd containers_build/boostdm/
$ sudo singularity build boostdm.sif Singularity
```
* Copy the singularity image in the root
```
$ cp containers_build/boostdm/boostdm.sif .
``` 

#### Nextflow

The pipeline runs with the [nextflow](https://www.nextflow.io/) workflow orchestrator. 
The boostDM-CH pipeline has been tested with nextflow version 20.07.1 which you can install e.g. with conda:
```
conda install -c bioconda nextflow=20.07.1
```

## Running the pipeline

#### Configuration

Consider the following adjustments before running the workflow:

* Adjust the values in the **nextflow.config** file. In particular, check the following:

   * env: **INTOGEN_DATASETS**, **DRIVERS_PATH** and **COHORTS_PATH** paths
   * singularity: **cacheDir** and **runOptions** paths.
   * process: make sure **container** matches the name of the singularity image file.
   * profile: include any custom execution profiles in the **config/** folder

* Adjust the **OUTPUT** parameter to match the **output/** folder in the **boostdm.nf** nextflow script.

* Execution profiles are defined in the **config/** folder. Make sure that your custom profile matches 
the computational resources available.

#### Run with nextflow

```
$ nextflow run boostdm.nf -profile <execution-profile> -resume
```

#### Run details

When the pipeline finishes, nextflow provides a summary of all the processes run in the file **trace.txt**. 
