# SARS-CoV-2_pipeline

Pipeline to map and analyze reads from SARS-CoV-2 NGS sequencing using Illumina Artic V3 protocol on Illumina NovaSeq at the H2030 Genome Center.

# Notes

This pipeline is designed to run within the H2030 Genome Center premises, but can be easily adapted to other infrastructures.

# Usage

./covpi.sh analysis.conf

## Important files

* __analysis.conf__: contains project and run specific information. The format is the following
```
export PROJECT=COVID19_160221
export RUN=210218_A00485_0104_AH2GNTDRXY
export RUNS_LANES=${RUN}:2
export OUT_DIR=/data/UHTS/2backup/projects/${PROJECT}/${RUN}
```

* __covpi.conf__: contains references in the local H2030 Genome Center environment. The values must be adapted to the environment 
```
###############################################################
# WORKING DIRS
###############################################################
export SRC_DIR=/data/UHTS/2backup/tools/covid_pipeline
export RUN_DIR=/scratch/permanent/PIPELINEOUTPUT/demultiplexing/

###############################################################
# REFERENCES
###############################################################
export REF_DIR=${SRC_DIR}
export SARS_REF=${REF_DIR}/refs/NC_045512/NC_045512.fasta
export REF=${REF_DIR}/refs/Homo_sapiens.GRCh38.99_NC_045512/Homo_sapiens.GRCh38.99_NC_045512.fa
export GTF=${REF_DIR}/refs/Homo_sapiens.GRCh38.99_NC_045512/Homo_sapiens.GRCh38.99_NC_045512.gtf
export INDEX=${REF_DIR}/refs/indexes/bowtie/Homo_sapiens.GRCh38.99_NC_045512
export VARIANT_ANNOTATION=${REF_DIR}/refs/NC_045512/variants_v1.csv

###############################################################
# BUILD INDEX
###############################################################
export BUILD_INDEX=0
```

__NB__: the reference directory containing the necessary files (references) for the execution of the pipeline is not distributed because of its size. 
You can ask the references to the to the contact person.

## Dependencies

The pipeline depends on common programs used for NGS data analysis:
* bowtie2
* samtools
* bedtools
* freeBayes
* R

## Contact

Please contact lorenzo.cerutti@health2030.ch or the H2030 Genome Center for questions and issues. 

