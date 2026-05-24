rule subsample_reads_for_assembly:
    input:
        r1="data/clean_fastq/{sample}_R1.fastq.gz",
        r2="data/clean_fastq/{sample}_R2.fastq.gz",
    output:
        r1=temp("data/assembly_reads/{sample}_R1.subsampled.fastq.gz"),
        r2=temp("data/assembly_reads/{sample}_R2.subsampled.fastq.gz"),
    log:
        "results/logs/assembly/{sample}.rasusa.log",
    threads:
        2
    params:
        genome_size=config["assembly"].get("estimated_genome_size", "6m"),
        coverage=config["assembly"].get("subsample_coverage", 120),
        seed=config["assembly"].get("subsample_seed", 42),
    conda:
        "../envs/spades_assembly.yaml"
    shell:
        """
        mkdir -p data/assembly_reads results/logs/assembly
        rasusa reads \
          {input.r1} \
          {input.r2} \
          --output {output.r1} \
          --output {output.r2} \
          --genome-size {params.genome_size} \
          --coverage {params.coverage} \
          --seed {params.seed} > {log} 2>&1
        """


rule spades_assembly:
    input:
        r1="data/assembly_reads/{sample}_R1.subsampled.fastq.gz",
        r2="data/assembly_reads/{sample}_R2.subsampled.fastq.gz",
    output:
        assembly="results/assemblies/{sample}/{sample}.assembly.fasta",
        scaffolds="results/assemblies/{sample}/{sample}.scaffolds.fasta",
        contigs="results/assemblies/{sample}/{sample}.contigs.fasta",
        spades_log="results/assemblies/{sample}/spades.log",
        params="results/assemblies/{sample}/params.txt",
    log:
        "results/logs/assembly/{sample}.spades.log",
    threads:
        config["resources"].get("assembly_threads", 8)
    params:
        memory_gb=config["assembly"].get("spades_memory_gb", 20),
        min_contig_length=config["assembly"].get("min_contig_length", 500),
        extra=config["assembly"].get("extra_args", "--isolate --only-assembler"),
        keep_workdir=config["assembly"].get("keep_spades_workdir", False),
        workdir=lambda wildcards: f"results/assemblies/{wildcards.sample}/spades_work",
    conda:
        "../envs/spades_assembly.yaml"
    shell:
        r"""
        mkdir -p results/assemblies/{wildcards.sample} results/logs/assembly
        rm -rf {params.workdir}
        spades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {params.workdir} \
          -t {threads} \
          -m {params.memory_gb} \
          {params.extra} > {log} 2>&1
        cp {params.workdir}/scaffolds.fasta {output.scaffolds}
        cp {params.workdir}/contigs.fasta {output.contigs}
        cp {params.workdir}/spades.log {output.spades_log}
        cp {params.workdir}/params.txt {output.params}
        python scripts/filter_fasta_by_length.py \
          --input {output.scaffolds} \
          --output {output.assembly} \
          --min-length {params.min_contig_length}
        if [ "{params.keep_workdir}" != "True" ] && [ "{params.keep_workdir}" != "true" ]; then
          rm -rf {params.workdir}
        fi
        """
