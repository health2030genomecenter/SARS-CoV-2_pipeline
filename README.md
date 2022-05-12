# SARS-CoV-2_pipeline
Pipeline to map and analyse reads from SARS-CoV-2 NGS sequencing using the Illumina Artic V4 protocol on Illumina NovaSeq6000 at the H2030 Genome Centre.

# Notes

This pipeline is designed to run in the H2030 Genome Centre facility, but can be easily adapted to other facilities. Ask the contact person if you need support.

# Usage

./covpipe.sh analysis.conf

## Important files

* __analysis.conf__: contains project and run specific information. The format is the following
```
export VERSION=GC_SARS-CoV-2_pipeline_v0.11
export RUN=220428_A00485_0274_AH3V2TDRX2
export PROJECT=COVID19_250422
export RUNS_LANES=220428_A00485_0274_AH3V2TDRX2:2
export PRJ_DIR=/data/UHTS/projects/COVID19/COVID19_250422/220428_A00485_0274_AH3V2TDRX2
export OUT_DIR=/data/UHTS/2backup/projects/COVID19_250422/220428_A00485_0274_AH3V2TDRX2
export FASTQ_DIR=${OUT_DIR}/fastq
```

* __covpipe.conf__: contains references in the local H2030 Genome Center environment. The values must be adapted to your environment 
```
###############################################################
# WORKING DIRS
###############################################################
export SRC_DIR=/data/UHTS/2backup/tools/covid_pipeline
export RUN_DIR=/scratch/permanent/PIPELINEOUTPUT/demultiplexing/
source ${SRC_DIR}/VERSION

###############################################################
# REFERENCES
###############################################################
export REF_DIR=${SRC_DIR}
export SARS_REF=${REF_DIR}/refs/NC_045512/NC_045512.fasta
export REF=${REF_DIR}/refs/Homo_sapiens.GRCh38.99_NC_045512/Homo_sapiens.GRCh38.99_NC_045512.fa
export GTF=${REF_DIR}/refs/Homo_sapiens.GRCh38.99_NC_045512/Homo_sapiens.GRCh38.99_NC_045512.gtf
export INDEX=${REF_DIR}/refs/indexes/bowtie/Homo_sapiens.GRCh38.99_NC_045512
export PRIMERS_BED=${REF_DIR}/refs/NC_045512/nCoV-2019.V4.bed
export VARIANT_ANNOTATION=${REF_DIR}/refs/NC_045512/variants_v8.csv
export HUMAN_COORDINATES=${REF_DIR}/refs/NC_045512/NC_045512_human_controls.coord

###############################################################
# External modules
###############################################################
export BOWTIE=UHTS/Aligner/bowtie2/2.3.4.1
export SAMTOOLS=samtools/1.13
export BCFTOOLS=UHTS/Samtools/1.10
export BEDTOOLS=UHTS/Analysis/BEDTools/2.26.0
export FREEBAYES=UHTS/Analysis/freebayes/1.2.0
export R_PACKAGE=R/3.6.1

###############################################################
# OTHER CONFS
###############################################################
export THREADS=4
export MIN_CVG=10

###############################################################
# BUILD INDEX
###############################################################
export BUILD_INDEX=0
```

__NB__: The reference directory containing the necessary files (references) for the execution of the pipeline is not distributed due to its size. 
You can request the references from the contact person and it will be sent to you. 

## Dependencies

The pipeline depends on common programs used for NGS data analysis:
* bowtie2
* samtools
* bedtools
* freeBayes
* R

## Contact

Please contact lorenzo.cerutti@health2030.ch or the H2030 Genome Centre with any questions or problems. 
