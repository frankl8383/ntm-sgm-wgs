rule prepare_assembly_qc_genomes:
    input:
        assemblies=expand("results/assemblies/{sample}/{sample}.assembly.fasta", sample=SAMPLES),
    output:
        genome_dir=directory("results/assembly_qc/genomes"),
    shell:
        """
        rm -rf {output.genome_dir}
        mkdir -p {output.genome_dir}
        for fasta in {input.assemblies}; do
          sample=$(basename "$fasta" .assembly.fasta)
          ln -sf "../../assemblies/${{sample}}/${{sample}}.assembly.fasta" "{output.genome_dir}/${{sample}}.fasta"
        done
        """


rule run_checkm2:
    input:
        genome_dir="results/assembly_qc/genomes",
        database=config["checkm2"]["database"],
    output:
        report="results/assembly_qc/checkm2/quality_report.tsv",
    log:
        "results/logs/checkm2.log",
    threads:
        config["resources"].get("checkm2_threads", 1)
    params:
        outdir="results/assembly_qc/checkm2",
        lowmem="--lowmem" if config.get("checkm2", {}).get("lowmem", True) else "",
        extra=config.get("checkm2", {}).get("extra_args", ""),
    conda:
        "../envs/checkm2.yaml"
    shell:
        """
        python scripts/run_checkm2_macos_fork.py \
          --input {input.genome_dir} \
          --extension .fasta \
          --output-directory {params.outdir} \
          --database-path {input.database} \
          --threads {threads} \
          {params.lowmem} \
          --force \
          --remove-intermediates \
          {params.extra} > {log} 2>&1
        """


rule run_quast_combined:
    input:
        assemblies=expand("results/assemblies/{sample}/{sample}.assembly.fasta", sample=SAMPLES),
    output:
        report_tsv="results/assembly_qc/quast_combined/report.tsv",
        transposed_tsv="results/assembly_qc/quast_combined/transposed_report.tsv",
        report_pdf="results/assembly_qc/quast_combined/report.pdf",
        report_html="results/assembly_qc/quast_combined/report.html",
    log:
        "results/logs/quast_combined.log",
    threads:
        config["resources"].get("quast_threads", 1)
    params:
        executable=config.get("quast", {}).get("executable", "quast.py"),
        labels=",".join(SAMPLES),
        outdir="results/assembly_qc/quast_combined",
        min_contig_length=config.get("quast", {}).get("min_contig_length", 500),
        extra=config.get("quast", {}).get("extra_args", ""),
    shell:
        """
        rm -rf {params.outdir}
        {params.executable} {input.assemblies} \
          -l {params.labels} \
          -o {params.outdir} \
          -t {threads} \
          -m {params.min_contig_length} \
          {params.extra} > {log} 2>&1
        """


rule run_gunc:
    input:
        genome_dir="results/assembly_qc/genomes",
        database=config["gunc"]["database"],
    output:
        report="results/assembly_qc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",
    log:
        "results/logs/gunc.log",
    threads:
        config["resources"].get("gunc_threads", 1)
    params:
        outdir="results/assembly_qc/gunc",
        extra=config.get("gunc", {}).get("extra_args", ""),
    conda:
        "../envs/gunc.yaml"
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        gunc run \
          -d {input.genome_dir} \
          -e .fasta \
          -r {input.database} \
          -o {params.outdir} \
          --temp_dir {params.outdir}/tmp \
          -t {threads} \
          {params.extra} > {log} 2>&1
        """


rule collect_assembly_qc:
    input:
        samplesheet=config["paths"]["samplesheet"],
        assemblies=expand("results/assemblies/{sample}/{sample}.assembly.fasta", sample=SAMPLES),
        quast_report="results/assembly_qc/quast_combined/transposed_report.tsv",
        checkm2_report="results/assembly_qc/checkm2/quality_report.tsv",
        gunc_report="results/assembly_qc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",
    output:
        summary="results/tables/assembly_qc_summary.tsv",
    log:
        "results/logs/collect_assembly_qc.log",
    params:
        thresholds="config/qc_thresholds.yaml",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/collect_assembly_qc.py \
          --samplesheet {input.samplesheet} \
          --assembly-dir results/assemblies \
          --thresholds {params.thresholds} \
          --quast-report {input.quast_report} \
          --checkm2-report {input.checkm2_report} \
          --gunc-report {input.gunc_report} \
          --output {output.summary} > {log} 2>&1
        """


rule integrate_read_assembly_judgement:
    input:
        read_judgement="results/tables/read_level_taxonomy_initial_judgement.tsv",
        assembly_qc="results/tables/assembly_qc_summary.tsv",
    output:
        "results/tables/read_assembly_initial_genome_level_judgement.tsv",
    log:
        "results/logs/integrate_read_assembly_judgement.log",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/integrate_read_assembly_judgement.py \
          --read-judgement {input.read_judgement} \
          --assembly-qc {input.assembly_qc} \
          --output {output} > {log} 2>&1
        """
