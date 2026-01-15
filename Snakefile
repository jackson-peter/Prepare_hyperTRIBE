# STAR alignment workflow
import glob
import os
import re

# Configuration
configfile: "config/config.yaml"

# Load config variables
#-------------------------------------------------------------------------------------------
RESULTS_DIR = config["results_dir"]
FASTQ_DIR = config["fastq_dir"]

GENOME_FASTA = config["genome_fasta"]
TRANSCRIPTOME_FASTA = config["transcriptome_fasta"]
GENOME_GTF = config["genome_gtf"]

# STAR parameters
STAR_INDEX_DIR = config["star_index_dir"]
STAR_SJDBOVERHANG = config["star_sjdbOverhang"]
STAR_GENOMESAINDEXNBASES = config["star_genomeSAindexNbases"]

# SALMON parameters
SALMON_INDEX = config["salmon_index_dir"]
SALMON_QUANT_DIR = os.path.join(RESULTS_DIR, "quantification", "salmon")
SALMON_THREADS = config["salmon_threads"]
SJDB_FILE = config.get("sjdb_file", "")  # Optional, empty string if not provided

# Output directories


LOG_DIR = os.path.join(RESULTS_DIR, "logs")
ALIGN_DIR = os.path.join(RESULTS_DIR, "alignment", "star")

MPILEUP_DIR = os.path.join(RESULTS_DIR, "variant", "mpileup")
TMP_DIR = os.path.join(RESULTS_DIR, "tmp")



R1_SUFFIX = config["r1_suffix"]  # e.g., "_1P.fq.gz"
R2_SUFFIX = R1_SUFFIX.replace("1", "2")  # Automatically deduce R2 suffix
PERL_SCRIPT = "scripts/RNAeditR_mpileup2bases.pl"

#-------------------------------------------------------------------------------------------

# Find all R1 files and extract sample names
r1_files = glob.glob(f"{FASTQ_DIR}/*{R1_SUFFIX}")
SAMPLES = []
for f in r1_files:
    basename = os.path.basename(f)
    # Remove R1 suffix to get sample name
    sample = basename.replace(R1_SUFFIX, '')
    SAMPLES.append(sample)

print(f"Found {len(SAMPLES)} samples: {SAMPLES}")
print(f"R1 suffix: {R1_SUFFIX}")
print(f"R2 suffix: {R2_SUFFIX}")

# Rule all
rule all:
    input:
        expand(f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(f"{SALMON_QUANT_DIR}/{{sample}}_quant", sample=SAMPLES),
        GENOME_FASTA + ".fai",
        STAR_INDEX_DIR,
        SALMON_INDEX,
        f"{MPILEUP_DIR}/all_regions_done.txt"


# STAR genome index generation
rule star_index:
    input:
        fasta = GENOME_FASTA,
        gtf = GENOME_GTF
    output:
        directory(STAR_INDEX_DIR)
    conda:
        "envs/star.yaml"
    params:
        genomeSAindexNbases = STAR_GENOMESAINDEXNBASES,
        sjdbOverhang = STAR_SJDBOVERHANG,
        tmp_base = TMP_DIR 
    threads: 20
    log:
        f"{LOG_DIR}/star_index.log"
    shell:
        """
        mkdir -p {params.tmp_base} $(dirname {log})
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
    conda:
        "envs/star.yaml"
    params:
        outprefix = lambda wildcards: f"{ALIGN_DIR}/{wildcards.sample}_",
        tmpdir = lambda wildcards: f"{TMP_DIR}/{wildcards.sample}",
        sjdb_param = f"--sjdbFileChrStartEnd {SJDB_FILE}" if SJDB_FILE else "",

    threads: 8

    log:
        f"{LOG_DIR}/star_align/{{sample}}.log"
    shell:
        """
        mkdir -p "$(dirname {log})" "{ALIGN_DIR}" "{TMP_DIR}"

        # Remove per-sample tmp dir if it exists (STAR requirement)
        rm -rf "{params.tmpdir}"

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

        """

# Salmon index generation rule
rule salmon_index:
    input:
        fasta = TRANSCRIPTOME_FASTA
    output:
        directory(SALMON_INDEX)
    conda:
        "envs/salmon.yaml"
    threads: 16
    log:
        f"{LOG_DIR}/salmon_index.log"
    shell:
        """   
        mkdir -p $(dirname {log})
        salmon index \
            -t {input.fasta} \
            -i {output} \
            -p {threads} \
            &> {log}
        """
# Salmon quantification rule
rule salmon_quant:
    input:
        r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
        r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX,
        index = SALMON_INDEX
    output:
        quant_dir = directory(f"{SALMON_QUANT_DIR}/{{sample}}_quant")
    conda:
        "envs/salmon.yaml"
    params:
        index = SALMON_INDEX
    threads: SALMON_THREADS
    log:
        f"{LOG_DIR}/salmon/{{sample}}.log"
    shell:
        """    
        mkdir -p $(dirname {log}) 
        salmon quant \
            -i {params.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            -o {output.quant_dir} \
            &> {log}
        """

# Index BAM files
rule index_bam:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam"
    output:
        bai = f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        "envs/samtools_perl.yaml"
    log:
        f"{LOG_DIR}/samtools_index/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        samtools index {input.bam}
        """


# Extract chromosome names
rule get_chromosomes:
    input:
        bam = f"{ALIGN_DIR}/{SAMPLES[0]}_Aligned.sortedByCoord.out.bam",  # Just first sample
        bai = f"{ALIGN_DIR}/{SAMPLES[0]}_Aligned.sortedByCoord.out.bam.bai"
    output:
        f"{MPILEUP_DIR}/chromosomes.txt"
    conda:
        "envs/samtools_perl.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        samtools idxstats {input.bam} | cut -f 1 | grep -v '*' > {output}
        """

rule faidx:
    input:
        fasta = GENOME_FASTA
    output:
        fai = GENOME_FASTA + ".fai"
    conda:
        "envs/samtools_perl.yaml"
    log:
        f"{LOG_DIR}/faidx.log"
    shell:
        """
        samtools faidx {input.fasta} 2> {log}
        """

# Generate regions for parallel processing
rule generate_regions:
    input:
        fasta = GENOME_FASTA,
        fai = GENOME_FASTA + ".fai",
        chroms = f"{MPILEUP_DIR}/chromosomes.txt"
    output:
        f"{MPILEUP_DIR}/regions.txt"
    conda:
        "envs/samtools_perl.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        while read chr; do
            samtools faidx {input.fasta} $chr | \
            awk -v chr="$chr" 'NR==1{{next}}{{len+=length($0)}}END{{for(i=0;i<len;i+=1000000){{start=i+1;end=(i+1000000>len)?len:i+1000000;print chr":"start"-"end}}}}' 
        done < {input.chroms} > {output}
        """

# Create BAM list for mpileup
rule create_bam_list:
    input:
        bams = expand(f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bais = expand(f"{ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)
    output:
        f"{MPILEUP_DIR}/bam_list.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        ls {input.bams} > {output}
        """

# Read regions and create a checkpoint
checkpoint split_regions:
    input:
        f"{MPILEUP_DIR}/regions.txt"
    output:
        directory(f"{MPILEUP_DIR}/regions_split")
    shell:
        """
        mkdir -p {output}
        split -l 1 {input} {output}/region_
        """

# Run mpileup for each region
rule mpileup_rnaeditr:
    input:
        bam_list = rules.create_bam_list.output[0],
        fasta = GENOME_FASTA,
        region_file = f"{MPILEUP_DIR}/regions_split/region_{{region_id}}"
    output:
        f"{MPILEUP_DIR}/output_{{region_id}}.txt"
    conda:
        "envs/samtools_perl.yaml"
    log:
        mpileup = f"{LOG_DIR}/mpileup/mpileup_{{region_id}}.log",
        perl = f"{LOG_DIR}/mpileup/perl_{{region_id}}.log"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """      
        mkdir -p $(dirname {log.mpileup}) $(dirname {output})
        region=$(cat {input.region_file})
        
        samtools mpileup \
            --max-depth 9999 \
            -f {input.fasta} \
            -r $region \
            -b {input.bam_list} \
            2> {log.mpileup} | \
        perl {PERL_SCRIPT} > {output} 2> {log.perl}
        """
# Aggregate all mpileup results
def aggregate_mpileup(wildcards):
    checkpoint_output = checkpoints.split_regions.get(**wildcards).output[0]
    region_ids = glob_wildcards(os.path.join(checkpoint_output, "region_{region_id}")).region_id
    return expand(f"{MPILEUP_DIR}/output_{{region_id}}.txt", region_id=region_ids)

rule aggregate_all_mpileup:
    input:
        aggregate_mpileup
    output:
        f"{MPILEUP_DIR}/all_regions_done.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        ls {MPILEUP_DIR}/output_*.txt | wc -l >> {output}
        """

