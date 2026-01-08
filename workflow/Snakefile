# STAR genome index generation rule

rule star_index:
    input:
        fasta = "TAIR10/Sequence/WholeGenomeFasta/WholeGenomeFastaWithADAR/genome_with_ADAR.fa",
        gtf = "TAIR10/Annotation/Genes/genes_with_ADAR.gtf"
    output:
        directory("index_STAR_ADAR")
    params:
        genomeSAindexNbases = 12,
        sjdbOverhang = 100
    threads: 20
    conda:
        "envs/star.yaml"
    log:
        "logs/star_index.log"
    shell:
        """
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
