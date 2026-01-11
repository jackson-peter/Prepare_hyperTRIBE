# STAR alignment workflow
import glob
import os
import re

# Configuration
configfile: "config/config.yaml"

# Load config variables
GENOME_FASTA = config["genome_fasta"]
GENOME_GTF = config["genome_gtf"]
STAR_INDEX_DIR = config["star_index_dir"]
STAR_SJDBOVERHANG = config["star_sjdbOverhang"]
STAR_GENOMESAINDEXNBASES = config["star_genomeSAindexNbases"]
ALIGN_DIR = config["align_dir"]
FASTQ_DIR = config["fastq_dir"]
TMP_DIR = config["tmp_dir"]
SJDB_FILE = config.get("sjdb_file", "")  # Optional, empty string if not provided
R1_SUFFIX = config["r1_suffix"]  # e.g., "_1P.fq.gz"
R2_SUFFIX = R1_SUFFIX.replace("1", "2")  # Automatically deduce R2 suffix

# Find all R1 files and extract sample names
r1_files = glob.glob(f"{FASTQ_DIR}/*{R1_SUFFIX}")
SAMPLES = []
for f in r1_files:
    basename = os.path.basename(f)
    # Remove R1 suffix to get sample name
    sample = basename.replace(R1_SUFFIX, '')
    SAMPLES.append(sample)

print(f"Found {len(SAMPLES)} samples for alignment: {SAMPLES}")
print(f"R1 suffix: {R1_SUFFIX}")
print(f"R2 suffix: {R2_SUFFIX}")

# Rule all
rule all:
    input:
        expand(f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)

# STAR genome index generation
rule star_index:
    input:
        fasta = GENOME_FASTA,
        gtf = GENOME_GTF
    output:
        directory(STAR_INDEX_DIR)
    params:
        genomeSAindexNbases = STAR_GENOMESAINDEXNBASES,
        sjdbOverhang = STAR_SJDBOVERHANG
    threads: 20
    # conda:
    #     "envs/star.yaml"
    log:
        "logs/star_index.log"
    shell:
        """
        export PATH=$PATH:/projects/renlab/apps/STAR/STAR-2.7.11a/source
        mkdir -p {output}
        STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeSAindexNbases {params.genomeSAindexNbases} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            2> {log}
        """

# STAR alignment rule
rule star_align:
    input:
        r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
        r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX,
        index = STAR_INDEX_DIR
    output:
        bam = f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam",
        log_final = f"{ALIGN_DIR}/{{sample}}_Log.final.out",
        log_progress = f"{ALIGN_DIR}/{{sample}}_Log.progress.out",
        log_out = f"{ALIGN_DIR}/{{sample}}_Log.out",
        sj = f"{ALIGN_DIR}/{{sample}}_SJ.out.tab",
        chimeric = f"{ALIGN_DIR}/{{sample}}_Chimeric.out.junction"
    params:
        outprefix = lambda wildcards: f"{ALIGN_DIR}/{wildcards.sample}_",
        tmpdir = lambda wildcards: f"{TMP_DIR}/{wildcards.sample}",
        sjdb_param = f"--sjdbFileChrStartEnd {SJDB_FILE}" if SJDB_FILE else ""
    threads: 80
    # conda:
    #     "envs/star.yaml"
    log:
        "logs/star_align/{sample}.log"
    shell:
        """
        export PATH=$PATH:/projects/renlab/apps/STAR/STAR-2.7.11a/source
        # Set ulimit
        ulimit -n 4096
        mkdir -p {TMP_DIR}
        # Remove temporary directory if it exists
        if [ -d "{params.tmpdir}" ]; then
            rm -rf "{params.tmpdir}"
        fi
        
        # Create output directory
        mkdir -p {ALIGN_DIR}
        
        mkdir -p $(dirname {log})
        
        # Run STAR
        STAR --runThreadN {threads} \
             --outTmpDir {params.tmpdir} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.outprefix} \
             {params.sjdb_param} \
             --chimJunctionOverhangMin 15 \
             --chimOutType WithinBAM SoftClip Junctions \
             --chimScoreMin 15 \
             --chimMultimapNmax 10 \
             --chimOutJunctionFormat 1 \
             --chimSegmentMin 12 \
             --outSAMattributes NH HI AS nM NM MD \
             --outFilterMismatchNoverReadLmax 0.4 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             --outFilterMismatchNmax 999 \
             --outFilterScoreMinOverLread 0.5 \
             --outFilterMatchNminOverLread 0.5 \
             2> {log}
        
        # Clean up temporary directory
        if [ -d "{params.tmpdir}" ]; then
            rm -rf "{params.tmpdir}"
        fi
        """