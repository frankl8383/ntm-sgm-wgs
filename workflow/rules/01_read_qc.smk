rule fastp_pe:
    input:
        r1=lambda wildcards: FASTQ_R1[wildcards.sample],
        r2=lambda wildcards: FASTQ_R2[wildcards.sample],
    output:
        r1="data/clean_fastq/{sample}_R1.fastq.gz",
        r2="data/clean_fastq/{sample}_R2.fastq.gz",
        html="results/qc/fastp/{sample}.fastp.html",
        json="results/qc/fastp/{sample}.fastp.json",
    log:
        "results/logs/fastp/{sample}.log",
    threads:
        config["resources"].get("fastp_threads", 4)
    params:
        qualified_quality_phred=config["fastp"].get("qualified_quality_phred", 20),
        length_required=config["fastp"].get("length_required", 50),
        adapter_flag="--detect_adapter_for_pe" if config["fastp"].get("detect_adapter_for_pe", True) else "",
        extra=config["fastp"].get("extra_args", ""),
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        mkdir -p results/qc/fastp results/logs/fastp data/clean_fastq
        fastp \
          --in1 {input.r1} \
          --in2 {input.r2} \
          --out1 {output.r1} \
          --out2 {output.r2} \
          --html {output.html} \
          --json {output.json} \
          --thread {threads} \
          --qualified_quality_phred {params.qualified_quality_phred} \
          --length_required {params.length_required} \
          {params.adapter_flag} \
          {params.extra} > {log} 2>&1
        """


rule fastqc_raw:
    input:
        lambda wildcards: FASTQ_R1[wildcards.sample] if wildcards.read == "R1" else FASTQ_R2[wildcards.sample]
    output:
        html="results/qc/fastqc/raw/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc/raw/{sample}_{read}_fastqc.zip",
    log:
        "results/logs/fastqc/raw/{sample}_{read}.log",
    threads:
        config["resources"].get("fastqc_threads", 1)
    params:
        outdir="results/qc/fastqc/raw",
        extra=config["fastqc"].get("extra_args", ""),
    conda:
        "../envs/fastqc_multiqc.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs/fastqc/raw
        fastqc --threads {threads} --outdir {params.outdir} {params.extra} {input} > {log} 2>&1
        """


rule fastqc_clean:
    input:
        lambda wildcards: f"data/clean_fastq/{wildcards.sample}_{wildcards.read}.fastq.gz"
    output:
        html="results/qc/fastqc/clean/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc/clean/{sample}_{read}_fastqc.zip",
    log:
        "results/logs/fastqc/clean/{sample}_{read}.log",
    threads:
        config["resources"].get("fastqc_threads", 1)
    params:
        outdir="results/qc/fastqc/clean",
        extra=config["fastqc"].get("extra_args", ""),
    conda:
        "../envs/fastqc_multiqc.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs/fastqc/clean
        fastqc --threads {threads} --outdir {params.outdir} {params.extra} {input} > {log} 2>&1
        """


rule multiqc_read_qc:
    input:
        preflight="results/tables/preflight_samplesheet_check.tsv",
        fastp_json=expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        raw_fastqc=expand("results/qc/fastqc/raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READS),
        clean_fastqc=expand("results/qc/fastqc/clean/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READS),
    output:
        html="results/qc/multiqc/multiqc_read_qc.html",
    log:
        "results/logs/multiqc_read_qc.log",
    params:
        outdir="results/qc/multiqc",
        title=config["multiqc"].get("title", "NTM SGM WGS read QC"),
    conda:
        "../envs/fastqc_multiqc.yaml"
    shell:
        """
        mkdir -p {params.outdir} results/logs
        multiqc \
          --title "{params.title}" \
          --filename multiqc_read_qc.html \
          --outdir {params.outdir} \
          results/qc/fastp results/qc/fastqc > {log} 2>&1
        """


rule collect_fastp_stats:
    input:
        samplesheet=config["paths"]["samplesheet"],
        fastp_json=expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
    output:
        "results/tables/fastp_read_qc_summary.tsv",
    log:
        "results/logs/collect_fastp_stats.log",
    params:
        thresholds="config/qc_thresholds.yaml",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        mkdir -p results/tables results/logs
        python scripts/collect_fastp_stats.py \
          --samplesheet {input.samplesheet} \
          --fastp-dir results/qc/fastp \
          --thresholds {params.thresholds} \
          --output {output} > {log} 2>&1
        """
