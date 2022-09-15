#!/usr/bin/bash

SELF=$(basename $0)
VERSION="2.1d"
THREADS=24
TF_DEDUPE=false

echo "Welcome to Honzo's Paired-end RNA-Seq pipeline ${SELF} v${VERSION}"

# HONZO'S NOTES
# VERSION 1.X
# developed pipeline that deduplicates with dedupe.sh, extracts UMIs with UMI tools
# pipeline capable of using both STAR or BWA (preferred) as an aligner
# workflow to process BAM/SAM files post-alignment
# calling on postProcess python script at the end
# VERSION 2.X
# "d" indicates dockerized version (directories changed)
# dockerizable version only contains human genome

# =====================================[ GETTING USER ARGUMENTS ]=====================================
# initial check
if [[ $# -eq 0 ]]; then
    echo "ERROR: No options have been chosen."
    echo "Please use RNASeqTargeted.sh -h for usage and options."
    exit 0
fi

# getting command line arguments from user
while [[ $# -gt 0 ]]; do
key="$1"
    case $key in
        -p|--project)
        PROJECT_ID="$2"
        shift # past argument
        shift # past value
        ;;
        -g|--genome)
        GENOME="$2"
        shift # past argument
        shift # past value
        ;;
        -d|--dedupe)
        TF_DEDUPE=true
        shift # past argument
        shift # past value
        ;;
        -h|--help)
        HELP=true
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
    esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters - whats this

if [[ -n $1 ]]; then
    echo "ERROR: no arguments were entered."
    echo "Please use ${SELF} -h for usage and options."
    exit 0
fi

if [[ "$HELP" = true ]]; then
    header.sh
	echo "USAGE: ${SELF} [-h] [-p] [-r]"
	echo ""
	echo "MAIN PARAMETERS:"
	echo -e "   -h, --help\t\tShow this message"
	echo -e "   -p, --project\tName of the project"
    echo -e "   -g, --genome\t\tSelect reference genome (e.g. hg38)"
    echo -e "   -d, --dedupe\t\tDeduplicate reads with dedupe.sh"
	echo ""
	exit 0
fi

# ======================================[ VARIABLE DECLARATIONS ]=====================================

# Setting directories and executable paths
DIR_PIPELINE="/app/"
DIR_WORK="/fastq/"

FASTA_HUMAN="/data/INDEX_HUMAN/GRCh38.p13_genomic.fna"
FASTA_TEST="/data/INDEX_TEST/test.fasta"

# Path variables: cutadapt
EXE_DEDUPE="/opt/bbmap/dedupe.sh"
EXE_REFORMAT="/opt/bbmap/reformat.sh"
EXE_PICARD="/opt/picard/build/libs/picard.jar"
EXE_BWA="/opt/bwa/bwa"

EXE_SAMTOOLS="/opt/samtools-1.16.1/samtools"
EXE_BCFTOOLS="/opt/bcftools-1.16/bcftools"

EXE_UMITOOLS="/usr/local/bin/umi_tools"
EXE_CUTADAPT="/usr/local/bin/cutadapt"

EXE_GATK="/opt/gatk-4.2.6.1/gatk"

#TruSeq Adapter Information
ADAPTER_1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

if [[ "${GENOME}" == "hg38" ]]; then
    INDEX_FASTA=${FASTA_HUMAN}

elif [[ "${GENOME}" == "test" ]]; then
    INDEX_FASTA=${FASTA_TEST}
fi



# =====================================[ CHECKING PRE-REQUISITES ]====================================

# check if index file exists, if not, create bwa index for the reference
echo "[rv] Checking for for genome indices"

# outputs:
DIR_REFERENCE=$(dirname ${INDEX_FASTA})
DICT_FASTA=${INDEX_FASTA/.f*a/.dict}

cd ${DIR_REFERENCE}
# check if the index files are there (only checking two)
if [ ! -f ${INDEX_FASTA}.amb ] || [ ! -f ${INDEX_FASTA}.sa ]; then
    echo "[rv] Reference FASTA missing index files. Indexing genome now..."
    ${EXE_BWA} index ${INDEX_FASTA}
fi
if [ ! -f ${INDEX_FASTA}.fai ]; then
    echo "[rv] Reference FASTA index not found. Indexing FASTA now..."
    ${EXE_SAMTOOLS} faidx ${INDEX_FASTA}
fi
if [ ! -f ${DICT_FASTA} ]; then
    echo "[rv] Reference FASTA dictionary not found. Creating sequence dictionary now..."
    java -jar ${EXE_PICARD} CreateSequenceDictionary R=${INDEX_FASTA} O=${DICT_FASTA}
fi

# =========================================[ START WORKFLOW ]=========================================

cd ${DIR_WORK}
RSV_LOG=${DIR_WORK}/RSV_commands.log

echo -e "[$(date +'%Y/%m/%d') $(date +'%T')] STARTING $(basename $0) for RNA-seq variant calling\n" > ${RSV_LOG}
echo -e "Working directory: ${DIR_WORK}"
echo -e "Reference genome: ${GENOME}"

[[ -d ${DIR_WORK}/RSV_OUTPUT ]] || mkdir ${DIR_WORK}/RSV_OUTPUT
for SAMPLE in $(ls ${DIR_WORK}/*.fastq.gz | xargs -n 1 basename | cut -f1 -d "R" | uniq); do
    
    FASTQ_RAW_1=${DIR_WORK}/${SAMPLE}R1_001.fastq.gz
    FASTQ_RAW_2=${DIR_WORK}/${SAMPLE}R2_001.fastq.gz
    SPL_NAME=$(echo ${FASTQ_RAW_1} | xargs -n 1 basename | cut -d '_' -f1,2)

    echo "--------------------------------------------------------------------------------"
    echo -e "Starting analysis for sample: ${SPL_NAME}"
    echo -e "\nWORKING ON SAMPLE ${SPL_NAME}" >> ${RSV_LOG}
    echo -e "Files: ${FASTQ_RAW_1}\t${FASTQ_RAW_2}" >> ${RSV_LOG}

    DIR_SAMPLE=${DIR_WORK}/RSV_OUTPUT/${SPL_NAME}
    [[ -d ${DIR_SAMPLE} ]] || mkdir ${DIR_SAMPLE}
    SPL_LOG=${DIR_SAMPLE}/${SPL_NAME}.log
    SPL_ERR=${DIR_SAMPLE}/${SPL_NAME}.err
    SPL_STAT=${DIR_SAMPLE}/${SPL_NAME}_readsummary.txt
    echo -e "Read Type\t${SPL_NAME}" > ${SPL_STAT}

    # Get number of raw reads, assuming number of seqs in R1 = number of seqs in R2
    echo "[rv] Getting number of raw reads in sample"
    NREADS_RAW=$(zcat ${FASTQ_RAW_1} | wc -l)/4 | bc
    echo -e "Raw FASTQ (x2)\t$(echo $(zcat ${FASTQ_RAW_1} | wc -l)/4 | bc)" >> ${SPL_STAT}

    # ------------------------------------[ Adapter Trimming ]----------------------------------------
    echo "[rv] Removing adapters from PE FASTQs with cutadapt..."
    # outputs:
    FASTQ_TRIMMED_1=${DIR_SAMPLE}/${SPL_NAME}_R1_trimmed.fastq.gz
    FASTQ_TRIMMED_2=${DIR_SAMPLE}/${SPL_NAME}_R2_trimmed.fastq.gz

    echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Running cutadapt:" >> ${RSV_LOG}
    echo -e "[cmd] cutadapt -j ${THREADS} --trim-n --max-n 0 --nextseq-trim=30 -m 20 -a ${ADAPTER_1} -A ${ADAPTER_2} -o ${FASTQ_TRIMMED_1} -p ${FASTQ_TRIMMED_2} ${FASTQ_RAW_1} ${FASTQ_RAW_2} >> ${SPL_LOG} 2>> ${SPL_ERR}" >> ${RSV_LOG}
    cutadapt -j ${THREADS} --trim-n --max-n 0 --nextseq-trim=30 -m 20 -a ${ADAPTER_1} -A ${ADAPTER_2} -o ${FASTQ_TRIMMED_1} -p ${FASTQ_TRIMMED_2} ${FASTQ_RAW_1} ${FASTQ_RAW_2} >> ${SPL_LOG} 2>> ${SPL_ERR}

    # Get number of adapter-trimmed reads
    echo -e "Adapter-trimmed (x2)\t$(echo $(zcat ${FASTQ_TRIMMED_1} | wc -l)/4 | bc)" >> ${SPL_STAT}

    # ----------------------------[ Collapsing reads (if applicable) ]--------------------------------

    if [ ${TF_DEDUPE} = true ]; then
        echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Deduplicating reads with BBTools" >> ${RSV_LOG}
        # dedupe.sh in1=read1.fq in2=read2.fq out=deduped.fq outd=duplicates ac=f
        FASTQ_DEDUPED=${DIR_SAMPLE}/${SPL_NAME}_deduped.fastq.gz
        echo -e "[cmd] ${EXE_DEDUPE} in1=${FASTQ_TRIMMED_1} in2=${FASTQ_TRIMMED_2} out=${FASTQ_DEDUPED} ac=f" >> ${RSV_LOG}
        ${EXE_DEDUPE} in1=${FASTQ_TRIMMED_1} in2=${FASTQ_TRIMMED_2} out=${FASTQ_DEDUPED} ac=f

        # reformat.sh in=deduped.fq out1=deduped1.fq out2=deduped2.fq
        FASTQ_DEDUPED_1=${DIR_SAMPLE}/${SPL_NAME}_R1_deduped.fastq.gz
        FASTQ_DEDUPED_2=${DIR_SAMPLE}/${SPL_NAME}_R2_deduped.fastq.gz
        echo -e "[cmd] ${EXE_REFORMAT} in=${FASTQ_DEDUPED} out1=${FASTQ_DEDUPED_1} out2=${FASTQ_DEDUPED_2}" >> ${RSV_LOG}
        ${EXE_REFORMAT} in=${FASTQ_DEDUPED} out1=${FASTQ_DEDUPED_1} out2=${FASTQ_DEDUPED_2}

        # Get number of deduplicated reads
        echo -e "Deduplicated (x2)\t$(echo $(zcat ${FASTQ_DEDUPED_1} | wc -l)/4 | bc)" >> ${SPL_STAT}

        # reassign FASTQ_TRIMMED variable name
        FASTQ_TRIMMED_1=${FASTQ_DEDUPED_1}
        FASTQ_TRIMMED_2=${FASTQ_DEDUPED_2}
    fi

    # -----------------------------[ Preliminary processing of UMIs ]---------------------------------

    echo -e "\n[rv] Replacing UMI sequences into FASTQ definitions with UMI-tools..."
    FASTQ_UMIEXTRACTED_1=${DIR_SAMPLE}/${SPL_NAME}_R1_umiextracted.fastq.gz
    FASTQ_UMIEXTRACTED_2=${DIR_SAMPLE}/${SPL_NAME}_R2_umiextracted.fastq.gz

    echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Running UMI extraction:" >> ${RSV_LOG}
    
    # UMI extract command
    # umi_tools extract --extract-method=string --bc-pattern=[PATTERN] --bc-pattern2=[PATTERN] 
    # --read2-in=[FASTQIN] --read2-out=[FASTQOUT] -L extract.log [OPTIONS]
    echo -e "[cmd] ${EXE_UMITOOLS} extract --stdin=${FASTQ_TRIMMED_1} --bc-pattern=NNNNNN --log=${SPL_LOG} --read2-in=${FASTQ_TRIMMED_2} --read2-out=${FASTQ_UMIEXTRACTED_2} --stdout ${FASTQ_UMIEXTRACTED_1}" >> ${RSV_LOG}
    ${EXE_UMITOOLS} extract --stdin=${FASTQ_TRIMMED_1} --bc-pattern=NNNNNN --log=${SPL_LOG} --read2-in=${FASTQ_TRIMMED_2} --read2-out=${FASTQ_UMIEXTRACTED_2} --stdout ${FASTQ_UMIEXTRACTED_1}
    
    # we don't really need this stat because it's going to be the same as Deduplicated (x2) given that all reads are long enough
    #echo -e "UMI-extracted (x2)\t$(echo $(zcat ${FASTQ_UMIEXTRACTED_1} | wc -l)/4 | bc)" >> ${SPL_STAT}

    # -----------------------------------[ Alignment to Genome ]--------------------------------------

    # outputs
    SAM_ALL=${DIR_SAMPLE}/${SPL_NAME}_genome_all.sam # alignment sam with mapped and unmapped reads
    BAM_ALLSORTED=${DIR_SAMPLE}/${SPL_NAME}_genome_all_sorted.bam
    TXT_IDXSTATS=${DIR_SAMPLE}/${SPL_NAME}_genome_indexStats.txt
    TXT_FLAGSTATS=${DIR_SAMPLE}/${SPL_NAME}_genome_flagStats.txt
    FASTQ_UNMAPPED_1=${DIR_SAMPLE}/${SPL_NAME}_R1_genome_unmapped.fastq.gz # unmapped reads in fastq format
    FASTQ_UNMAPPED_2=${DIR_SAMPLE}/${SPL_NAME}_R2_genome_unmapped.fastq.gz # unmapped reads in fastq format
    BAM_ALIGNED=${DIR_SAMPLE}/${SPL_NAME}_genome_Aligned.bam

    echo "[rv] Aligning UMI-extracted reads to genome with BWA"
    echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Aligning genome to reference with BWA mem" >> ${RSV_LOG}

    # verbose level (-v) 2 for writing all warnings and messages
    echo -e "[cmd] ${EXE_BWA} mem -t ${THREADS} -v 2 ${INDEX_FASTA} ${FASTQ_UMIEXTRACTED_1} ${FASTQ_UMIEXTRACTED_2} > ${SAM_ALL} 2>> ${SPL_ERR}" >> ${RSV_LOG}
    ${EXE_BWA} mem -t ${THREADS} -v 2 ${INDEX_FASTA} ${FASTQ_UMIEXTRACTED_1} ${FASTQ_UMIEXTRACTED_2}  > ${SAM_ALL} 2>> ${SPL_ERR}

    # getting mapping statistics
    echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Getting mapping summary and statistics" >> ${RSV_LOG}
    echo -e "[cmd] ${EXE_SAMTOOLS} sort -@ ${THREADS} -o ${BAM_ALLSORTED} ${SAM_ALL}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} sort -@ ${THREADS} -o ${BAM_ALLSORTED} ${SAM_ALL}
    echo -e "[cmd] ${EXE_SAMTOOLS} index ${BAM_ALLSORTED}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} index ${BAM_ALLSORTED}
    echo -e "[cmd] ${EXE_SAMTOOLS} idxstats ${BAM_ALLSORTED} > ${TXT_IDXSTATS}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} idxstats ${BAM_ALLSORTED} > ${TXT_IDXSTATS}
    echo -e "[cmd] ${EXE_SAMTOOLS} flagstat -@ ${THREADS} ${BAM_ALLSORTED} > ${TXT_FLAGSTATS}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} flagstat -@ ${THREADS} ${BAM_ALLSORTED} > ${TXT_FLAGSTATS}
    #rm ${SAM_ALL}

    # http://broadinstitute.github.io/picard/explain-flags.html 
    # we're keeping only mapped reads (doesn't matter if mate is unmapped)
    echo -e "[cmd] ${EXE_SAMTOOLS} view -b -F4 ${BAM_ALLSORTED} > ${BAM_ALIGNED}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} view -b -F4 ${BAM_ALLSORTED} > ${BAM_ALIGNED}

    # reads that we want to filter out (unmapped and improper-mapping)
    echo -e "${EXE_SAMTOOLS} fastq -F67 ${BAM_ALLSORTED} | gzip > ${FASTQ_UNMAPPED_1}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} fastq -F67 ${BAM_ALLSORTED} | gzip > ${FASTQ_UNMAPPED_1}
    echo -e "${EXE_SAMTOOLS} fastq -F131 ${BAM_ALLSORTED} | gzip > ${FASTQ_UNMAPPED_2}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} fastq -F131 ${BAM_ALLSORTED} | gzip > ${FASTQ_UNMAPPED_2}

    # get number of mapped reads
    echo -e "All mapped\t$(${EXE_SAMTOOLS} view -c ${BAM_ALIGNED})" >> ${SPL_STAT}
    echo -e "Proper-paired mapped\t$(${EXE_SAMTOOLS} view -c -f3 ${SAM_ALL})" >> ${SPL_STAT}


    # ----------------------------------[ Processing Alignments ]-------------------------------------
    #outputs:
    BAM_VALIDATED=${DIR_SAMPLE}/${SPL_NAME}_genome_aligned_validated.bam

    echo "[rv] Processing Alignment files..."
    echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Processing BAM files with samtools" >> ${RSV_LOG} 

    # echo "[rv] Sorting alignment file: ${BAM_ALIGNED}"
    # echo -e "[cmd] ${EXE_SAMTOOLS} sort -@ ${THREADS} -o ${BAM_SORTED} ${BAM_ALIGNED}" >> ${RSV_LOG}
    # ${EXE_SAMTOOLS} sort -@ ${THREADS} -o ${BAM_SORTED} ${BAM_ALIGNED}

    # Validate the BAM file
    echo "[rv] Validating BAM file..."
    echo -e "[CMD] java -jar ${EXE_PICARD} ValidateSamFile -M SUMMARY -IGNORE_WARNINGS true -I ${BAM_ALIGNED} -R ${INDEX_FASTA}" >> ${RSV_LOG}
    java -jar ${EXE_PICARD} ValidateSamFile -M SUMMARY -IGNORE_WARNINGS true -I ${BAM_ALIGNED} -R ${INDEX_FASTA} >> ${SPL_LOG}
    
    if grep -q ERROR:MISSING_READ_GROUP ${SPL_LOG}; then
        echo "Adding read group to BAM file for use with GATK"
        echo "[cmd] ${EXE_GATK} AddOrReplaceReadGroups -I ${BAM_ALIGNED} -O ${BAM_VALIDATED} -SORT_ORDER coordinate -LB lib1 -PL ILLUMINA -PU 1 -SM ${SPL_NAME}" >> ${RSV_LOG}
        ${EXE_GATK} AddOrReplaceReadGroups -I ${BAM_ALIGNED} -O ${BAM_VALIDATED} -SORT_ORDER coordinate -LB lib1 -PL ILLUMINA -PU 1 -SM ${SPL_NAME}

    elif ! grep "No errors found" ${SPL_LOG}; then
        echo "[Warning] Found unknown error in sample: ${SPL_NAME}"
        echo -e "[Warning] Unable to validate file: ${BAM_ALIGNED}\n\tMoving on to next sample..." >> ${RSV_LOG}
        continue

    else
        echo "[rv] Validation complete. No errors found."
        mv ${BAM_ALIGNED} ${BAM_VALIDATED}
    fi

    echo -e "[cmd] ${EXE_SAMTOOLS} index ${BAM_SORTED}" >> ${RSV_LOG}
    ${EXE_SAMTOOLS} index ${BAM_VALIDATED}

    # ------------------------------------[ Calling Variants ]----------------------------------------

    #outputs:
    VCF_FILE=${DIR_SAMPLE}/${SPL_NAME}_genome_variants.vcf
    TXT_VCFSTATS=${DIR_SAMPLE}/${SPL_NAME}_genome_variantStats.txt

    echo "[rv] Calling variants with HaplotypeCaller from GATK..."
    echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Calling variants with HaplotypeCaller" >> ${RSV_LOG}

    echo -e "[cmd] ${EXE_GATK} --java-options "-Xmx4g" HaplotypeCaller -R ${INDEX_FASTA} -I ${BAM_VALIDATED} -O ${VCF_FILE}" >> ${RSV_LOG}
    ${EXE_GATK} --java-options "-Xmx4g" HaplotypeCaller -R ${INDEX_FASTA} -I ${BAM_VALIDATED} -O ${VCF_FILE}

    echo -e "[cmd] ${EXE_BCFTOOLS} stats ${VCF_FILE} > ${TXT_VCFSTATS}" >> ${RSV_LOG}
    ${EXE_BCFTOOLS} stats ${VCF_FILE} > ${TXT_VCFSTATS}

done

# ==========================================[ POST PROCESS ]==========================================

echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Running post-process" >> ${RSV_LOG}
echo -e "[cmd] python3 ${DIR_PIPELINE}/rsv_postprocess.py ${PROJECT_ID}" >> ${RSV_LOG}
python3 ${DIR_PIPELINE}/rsv_postprocess.py ${PROJECT_ID}

echo -e "\n[$(date +'%Y/%m/%d') $(date +'%T')] Program complete." >> ${RSV_LOG}
echo "[rv] Program complete."