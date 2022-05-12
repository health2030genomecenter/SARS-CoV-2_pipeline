#!/bin/env bash

PARTITION=normal #queue10A|queue10B|normal
EXEC_DIR=/data/UHTS/2backup/tools/SARS-CoV-2_pipeline/curr
source ${EXEC_DIR}/covpipe.conf
source ${SRC_DIR}/src/utils.sh
uid=$(uid_gen)

# READ INPUT CONFIG FILE passed as argument
source $1 
if [ "x${PROJECT}" == "x" ]; then
    echo "Missing analysis configuration file"
    exit 1
fi


###############################################################
# SETUP WORKING DIR
###############################################################
mkdir -p ${PRJ_DIR}
ln -s ${PRJ_DIR} ${OUT_DIR}
cd ${OUT_DIR}
echo $VERSION > VERSION

###############################################################
# BUILD BOWTIE2 INDEX IF REQUESTED
#############################################################
if [ "x1" == "x$BUILD_INDEX" ]; then
    mkdir -p $INDEX
    module add $BOWTIE
    bowtie2-build $REF $INDEX
    module remove $BOWTIE
fi


###############################################################
# GET RUNS
###############################################################
runs_lanes=$(echo ${RUNS_LANES} | tr "," "\n")


###############################################################
# COLLECT FASTQs
###############################################################
echo_info "PART0 - Collect fastq files and associate to sample"
declare -A samples
declare -A samples2
for run_lane in ${runs_lanes}; do
    runName=$(echo ${run_lane} | cut -d ':' -f 1)
    totLanes=$(echo ${run_lane} | cut -d ':' -f 2)
    run=${RUN_DIR}/${runName}
    for lane in `seq 1 $totLanes`; do
        for i in `ls ${run}/demult1/part$((lane-1))/${PROJECT}/*_L00${lane}_R1_001.fastq.gz 2>/dev/null`; do 
	        a=$(echo $i | sed "s/_L00${lane}_R1_001.fastq.gz//" | xargs basename | sed "s/_null_Capt.*//")
            # Collect R1s
            for j in `ls ${run}/demult1/part$((lane-1))/${PROJECT}/${a}_null_*_L00${lane}_R1_*.fastq.gz`; do
                if [ "x"${samples[$(to_string $a)]} != "x" ]; then
                    samples[$(to_string $a)]=${samples[$(to_string $a)]}",$j" # concatenate file names
                else
                    samples[$(to_string $a)]=$j # first file name
                fi
            done
            # Collect R2s
            for j in `ls ${run}/demult1/part$((lane-1))/${PROJECT}/${a}_null_*_L00${lane}_R2_*.fastq.gz`; do
                if [ "x"${samples2[$(to_string $a)]} != "x" ]; then
                    samples2[$(to_string $a)]=${samples2[$(to_string $a)]}",$j"
                else
                    samples2[$(to_string $a)]=$j
                fi
            done
        done
    done
done


###############################################################
# MAPPING
###############################################################
echo_info "PART1 - Mapping of reads to reference"
samdir=sam
jobname=part1_mapping.$uid
if [ "x1" == "x"$(pi_lock $jobname) ]; then
    set_dir ${samdir}
    for sample in "${!samples[@]}"; do
        out=${samdir}/${sample}.sam
        jobid=$(sbatch --parsable -p ${PARTITION} --mem 80G --nodes 1 -n 20 -o sam_logs/${sample}.out -e sam_logs/${sample}.err --job-name $jobname --wrap="\ 
            module add $BOWTIE; \
                bowtie2 -x ${INDEX} -1 ${samples[$(to_string $sample)]} -2 ${samples2[$(to_string $sample)]} -S $out -p ${THREADS} --no-mixed --no-unal --no-discordant --sensitive --rdg 5,1.9")
        echo "... processing $sample JID=$jobid"
        echo "-1 ${samples[$(to_string $sample)]}"
        echo "-2 ${samples2[$(to_string $sample)]}"
    done
fi
pi_wait $jobname

###############################################################
# BAM SORTING AND INDEXING
############################################################### 
echo_info "PART2 - Sort and index BAMs"
jobname=part2_bamSorting.$uid
bamdir=bam
if [ "x1" == "x"$(pi_lock $jobname) ]; then
    set_dir bam
    for sample in "${!samples[@]}"; do
        s="${samdir}/${sample}.sam"
        b="${bamdir}/${sample}.bam"
        c="${bamdir}/${sample}"
        jobid=$(sbatch --parsable -p ${PARTITION} --mem 80G --nodes 1 -n 20 -o bam_logs/${sample}.out -e bam_logs/${sample}.err --job-name $jobname --wrap="\
            module add $SAMTOOLS; \
            samtools ampliconclip -@ ${THREADS} --hard-clip -u -b ${PRIMERS_BED} -o ${c}.cl.bam ${s}; \ 
            samtools fixmate -@ ${THREADS} -r -m ${c}.cl.bam ${c}.fm.bam; \ 
            samtools sort ${c}.fm.bam -@ ${THREADS} -O BAM -o ${b}; \
            samtools index -@ ${THREADS} ${b} > ${b}.bai")
        echo "... processing $sample JID=$jobid"
    done 
fi
pi_wait $jobname;

###############################################################
# CALCULATE  COVERAGE
############################################################### 
echo_info "PART3 - Extract coverage"
jobname=part3_coverage.$uid
if [ "x1" == "x"$(pi_lock $jobname) ]; then
    set_dir coverage_expanded
    for sample in "${!samples[@]}"; do
	    b="${bamdir}/${sample}.bam"
	    k="coverage_expanded/${sample}.cvg"
        jobid=$(sbatch --parsable -p ${PARTITION} --mem 80G --nodes 1 -n 20 -o coverage_expanded_logs/${sample}.out -e coverage_expanded_logs/${sample}.err --job-name $jobname --wrap="\
		    module add $SAMTOOLS; \
		    samtools depth -J -a -r NC_045512:1-29903 $b > $k")
        echo "... processing $sample JID=$jobid"
    done	
fi
pi_wait $jobname


###############################################################
# CALCULATE COVERAGE SUMMARY
############################################################### 
echo_info "PART4 - Get coverage summary"
jobname=part4_coverageSummary.$uid
if [ "x1" == "x"$(pi_lock $jobname) ]; then
    set_dir coverage_summary
    for sample in "${!samples[@]}"; do
	    b="${bamdir}/${sample}.bam"
	    k="coverage_summary/${sample}.txt"
	    jobid=$(sbatch --parsable -p ${PARTITION} --mem 80G --nodes 1 -n 20 -o coverage_summary_logs/${sample}.out -e coverage_summary_logs/${sample}.err --job-name $jobname --wrap="\
		    module add $SAMTOOLS; \
		    samtools coverage -r NC_045512:1-29903 $b > $k")
        echo "... processing $k JID=$jobid"
    done
fi
pi_wait $jobname # wait COVERAGE SUMMARY step4

###############################################################
# CALCULATE COVERAGE EXPANDED REMOVED (Part5) AS DONE IN Part3
############################################################### 
echo_info "PART5 - Expand coverage (REMOVED)"

###############################################################
#  BUILD VCF (start it now as it is time consuming)
############################################################### 
echo_info "PART6 - Build VCFs"
jobname=part6_variantMap.$uid
if [ "x1" == "x"$(pi_lock $jobname) ]; then
    set_dir vcf
    for sample in "${!samples[@]}"; do
        b="${bamdir}/${sample}.bam"
        v="vcf/${sample}.vcf"
        jobid=$(sbatch --parsable -p ${PARTITION} --mem 80G --nodes 1 -n 20 -o vcf_logs/${sample}.out -e vcf_logs/${sample}.err --job-name $jobname --wrap="\
	        module add $FREEBAYES; \
	        freebayes -p 1 -j --region NC_045512:1..29903 -f $REF $b > $v;")
        echo "... processing $x JID=$jobid"
    done
fi
#pi_wait $jobname # postpone as we can proceed with other tasks


###############################################################
# COUNTS AND TABLES
###############################################################
echo_info "PART7 - Coverage report"
localname=part7_coverageReport.$uid
if [ "x1" == "x"$(pi_lock $localname) ]; then
    echo_info "PART7a - Counts from FASTQs"
    for run_lane in ${runs_lanes}; do
        runName=$(echo ${run_lane} | cut -d ':' -f 1)
        stat_json="${stat_json} ${RUN_DIR}/${runName}/demult1/Stats.json"
    done
    $SRC_DIR/src/get_tot_reads.pl $stat_json > counts.csv

    echo_info "PART7b - Counts from BAMs"
    module add $SAMTOOLS
    cat counts.csv | while read i; do
	    a=$(echo $i | cut -d ',' -f 1 | sed 's/_null_Capt.*$//')
	    n=$(echo $i | cut -d ',' -f 2)
        if [[ "x" != "x"${samples[$(to_string $a)]} ]]; then
            if [ -f ${bamdir}/${a}.bam ]; then 
		        h=0
                while read j; do
                    t=0
                    t=$(samtools stats bam/${a}.bam ${j} | grep 'reads mapped:' | awk '{print $4}')
                    h=$((h+t))
                done <<< $(cat ${HUMAN_COORDINATES})              
                
                c=$(samtools view -@ 8 -c ${bamdir}/${a}.bam)
		        echo "$a,$n,$c,$h"
            else
                echo >&2 "ERROR missing ${bamdir}/${a}.bam file"
                exit 2
            fi
        fi
    done > counts2.csv
    module remove $SAMTOOLS

    echo_info "PART7c - Build summary table"
    echo 'sample,tot_numreads,tot_mapreads,human_bait_reads,refname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq' > summaryR.csv
    cat counts2.csv | while read l; do 
	    n=$(echo $l | cut -d ',' -f 1)
	    c=$(tail -1 coverage_summary/${n}.txt)
	    echo "$l,$c"
    done | sed 's/\t/,/g' >> summaryR.csv

    echo_info "PART7d - Build summary table with >= 15x"
    export CVG=15
    echo 'sample,tot_numreads,tot_mapreads,human_bait_reads,refname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq,numpos15x,perpos15x' > summaryR2.csv
    cat summaryR.csv | grep -v 'sample' | while read l; do 
        name=$(echo $l | cut -d ',' -f 1)
	    n=$(cat coverage_expanded/${name}.cvg | awk -v x=$CVG '{if ($3 >= x) {print}}' | wc -l)
	    p=$(echo "scale=2;100*$n/29903" | bc)
	    echo "$l,$n,$p"
    done >> summaryR2.csv

    echo_info "PART7e - Build summary table with >= 50x"
    export CVG=50
    echo 'sample,tot_numreads,tot_mapreads,human_bait_reads,refname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq,numpos15x,perpos15x,numpos50x,perpos50x' > summaryR3.csv
    cat summaryR2.csv | grep -v 'sample' | while read l; do 
	    name=$(echo $l | cut -d ',' -f 1)
	    n=$(cat coverage_expanded/${name}.cvg |  awk -v x=$CVG '{if ($3 >= x) {print}}' | wc -l)
	    p=$(echo "scale=2;100*$n/29903" | bc)
	    echo "$l,$n,$p"
    done >> summaryR3.csv

    echo_info "PART7f - Build summary table with < 5x"
    export CVG=5
    echo 'sample,tot_numreads,tot_mapreads,human_bait_reads,refname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq,numpos15x,perpos15x,numpos50x,perpos50x,perN5' > summaryR4.csv
    cat summaryR3.csv | grep -v 'sample' | while read l; do
	    name=$(echo $l | cut -d ',' -f 1)
	    n=$(cat coverage_expanded/${name}.cvg | awk -v x=$CVG '{if ($3 < x) {print}}' | wc -l)
	    p=$(echo "scale=2;100*$n/29903" | bc)
	    echo "$l,$p"
    done >> summaryR4.csv

    #echo_info "PART7g - Build consensi mask for < 5x"
    #for j in `ls coverage_expanded/`; do 
    #    k=${j/\.cvg/.mask}
    #    cat coverage_expanded/$j | awk '{if ($2 <= 5) { print "NC_045512\t"$1+1 }}' > vcf/$k
    #done

    echo_info "PART7h - Coverage report"
    d=$(date)
    echo "---" > R.Rmd
    echo "title: \"$PROJECT\"" >> R.Rmd
    echo "date: \"$d\"" >> R.Rmd
    cat ${SRC_DIR}/src/R.Rmd >> R.Rmd 
    #module add R/latest
    module add $R_PACKAGE
    Rscript -e "rmarkdown::render('R.Rmd')"
    module remove $R_PACKAGE 
    mv -v R.html ${PROJECT}.html
fi
pi_unlock $localname
   

###############################################################
# Wait for VCF generation terminates
###############################################################
pi_wait $jobname # wait VCF build


###############################################################
# Build consensi
###############################################################
echo_info "PART8 - Filter VCFs and build consensi"
jobname=part8_buildConsensi.$uid
if [ "x1" == "x"$(pi_lock $jobname) ]; then
    set_dir consensi
    for sample in "${!samples[@]}"; do
        v="vcf/${sample}.vcf"
        c="coverage_expanded/${sample}.cvg"
        f="consensi/tmp_${sample}.fasta"
        z="vcf/${sample}.vcf.gz"
        m="vcf/${sample}.mask"
        jobid=$(sbatch --parsable -p ${PARTITION} --mem 80G --nodes 1 -n 20 -o consensi_logs/${sample}.out -e consensi_logs/${sample}.err --job-name $jobname --wrap="\
	        module add $BCFTOOLS; \
	        bcftools norm --threads ${THREADS} -m -any $v | 
            bcftools view --threads ${THREADS} -i 'QUAL > 650 || (QUAL > 1 && QUAL / INFO/AO > 6 && (GL[-:1] - GL[-:0]) > 10 )' -Oz -o $z; \
	        bcftools index --threads ${THREADS} $z; \
	        zcat ${z} | $SRC_DIR/src/buildConsensus.pl -m $MIN_CVG -f ${SARS_REF} -c ${c} > ${f};")
    done
fi
pi_wait $jobname
            
#bcftools view --threads ${THREADS} -i 'QUAL > 650 || (QUAL > 1 && QUAL / INFO/AO > 10 && DP > 5 && SAF > 0 && SAR > 0 && (GL[-:1] - GL[-:0]) > 10 )' -Oz -o $z; \


###############################################################
# Post-processing of consensi
###############################################################
echo_info "PART9 - Consensi post-processing"
localname=part9_consensiPostprocessing.$uid
if [ "x1" == "x"$(pi_lock $localname) ]; then
    echo_info "PART9a - Rename the consensi"
    for sample in "${!samples[@]}"; do
	    f=${sample}.fasta; \
	    cat consensi/tmp_${f} | sed 's/^>/>'"${sample}"' /' > consensi/${f}; \
    done

    echo_info "PART9b - Select consensi where 95% > 15x"
    set_dir consensi_gisaid_ok
    cat summaryR4.csv  | awk -F ',' '{if ($15 >= 95) { print $1 }}'  | grep -v sample | while read l; do cp -v consensi/${l}.fasta consensi_gisaid_ok/; done

    echo_info "PART9c - Select consensi where 95% > 15x is not reached"
    set_dir consensi_gisaid_not_ok
    cat summaryR4.csv  | awk -F ',' '{if ($15 < 95) { print $1 }}'  | grep -v sample | while read l; do cp -v consensi/${l}.fasta consensi_gisaid_not_ok/; done
fi
pi_unlock $localname


###############################################################
# Post-processing of variants 
###############################################################
echo_info "PART10 - Variants post-processing"
localname=part10_variantsPostprocessing.$uid
if [ "x1" == "x"$(pi_lock $localname) ]; then
    echo_info "PART10a - Collect mutations"
    module add $BCFTOOLS
    for sample in "${!samples[@]}"; do
	    bcftools view vcf/${sample}.vcf.gz | ${SRC_DIR}/src/vcf2tbl.pl -n $sample
    done > vcf.csv
    module remove $SAMTOOLS

    echo_info "PART10b - Variant of concern attribution"
    cat vcf.csv | ${SRC_DIR}/src/annvar.pl -m $MIN_CVG -c 0.6 -d coverage_expanded/ -v $VARIANT_ANNOTATION  2>variant_attribution.log > variant_attribution.csv
fi
pi_unlock $localname
