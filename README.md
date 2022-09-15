# RSV
Basic Variant-Calling pipeline with Read-mapping Diagnostics

## HONZO'S NOTES
### VERSION 1.X
- created flag to deduplicate reads with with dedupe.sh, extracts UMIs with UMI tools
- pipeline capable of using both STAR or BWA (preferred) as an aligner
- workflow to process BAM/SAM files post-alignment
- calling on postProcess python script at the end
### VERSION 2.X
- "d" indicates dockerized version (directories changed)
- small updates to postProcess.py script
### FUTURE UPDATES
- ADAPTER_1 and ADAPTER_2 sequences are currently hardcoded, should be reading from file instead
- Currently two genomes set up, hg38 and angle, planning to add more
- We should look into what kind of plots we can generate for postProcess
- Perform FastQC on reads (already in docker image)

## DEPENDENCIES
- Picard tools (https://github.com/broadinstitute/picard.git)
- GATK (available locally as zipped file)
- FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip)
- Samtools, Bcftools (available locally as zipped file)
- BWA (https://github.com/lh3/bwa.git)
- Python libraries listed in requirements.txt
