rule kraken2_clean:
    input:
        r1="data/clean_fastq/{sample}_R1.fastq.gz",
        r2="data/clean_fastq/{sample}_R2.fastq.gz",
    output:
        report="results/taxonomy/kraken2/{sample}.kraken2.report.txt",
        classifications="results/taxonomy/kraken2/{sample}.kraken2.out.gz",
    log:
        "results/logs/kraken2/{sample}.log",
    threads:
        config["resources"].get("kraken2_threads", 8)
    params:
        db=config["kraken2"].get("database"),
        confidence=config["kraken2"].get("confidence", 0.05),
        minimum_base_quality=config["kraken2"].get("minimum_base_quality", 0),
        memory_mapping="--memory-mapping" if config["kraken2"].get("memory_mapping", True) else "",
        save_per_read_output=config["kraken2"].get("save_per_read_output", False),
        extra=config["kraken2"].get("extra_args", ""),
    conda:
        "../envs/kraken_bracken.yaml"
    shell:
        r"""
        mkdir -p results/taxonomy/kraken2 results/logs/kraken2
        if [ -z "{params.db}" ] || [ "{params.db}" = "None" ] || [ "{params.db}" = "null" ]; then
          echo "Kraken2 database is not configured. Set kraken2.database in config/config.yaml." > {log}
          exit 1
        fi
        if [ "{params.save_per_read_output}" = "True" ] || [ "{params.save_per_read_output}" = "true" ]; then
          kraken2 \
            --db {params.db} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --use-names \
            --confidence {params.confidence} \
            --minimum-base-quality {params.minimum_base_quality} \
            {params.memory_mapping} \
            {params.extra} \
            --report {output.report} \
            --output - \
            {input.r1} {input.r2} 2> {log} | gzip -c > {output.classifications}
        else
          kraken2 \
            --db {params.db} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --use-names \
            --confidence {params.confidence} \
            --minimum-base-quality {params.minimum_base_quality} \
            {params.memory_mapping} \
            {params.extra} \
            --report {output.report} \
            --output /dev/null \
            {input.r1} {input.r2} > /dev/null 2> {log}
          printf "" | gzip -c > {output.classifications}
        fi
        """


rule bracken_species:
    input:
        report="results/taxonomy/kraken2/{sample}.kraken2.report.txt",
    output:
        abundance="results/taxonomy/bracken/{sample}.species.bracken.tsv",
        report="results/taxonomy/bracken/{sample}.species.kraken2.report.txt",
    log:
        "results/logs/bracken/{sample}.species.log",
    threads:
        config["resources"].get("bracken_threads", 1)
    params:
        db=config["bracken"].get("database") or config["kraken2"].get("database"),
        read_length=config["bracken"].get("read_length", 150),
        threshold=config["bracken"].get("threshold", 10),
        level=config["bracken"].get("species_level", "S"),
        extra=config["bracken"].get("extra_args", ""),
    conda:
        "../envs/kraken_bracken.yaml"
    shell:
        """
        mkdir -p results/taxonomy/bracken results/logs/bracken
        if [ -z "{params.db}" ] || [ "{params.db}" = "None" ] || [ "{params.db}" = "null" ]; then
          echo "Bracken database is not configured. Set bracken.database in config/config.yaml." > {log}
          exit 1
        fi
        if command -v bracken >/dev/null 2>&1; then
          bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.abundance} \
            -w {output.report} \
            -r {params.read_length} \
            -l {params.level} \
            -t {params.threshold} \
            {params.extra} > {log} 2>&1
        elif command -v est_abundance.py >/dev/null 2>&1; then
          kmer_distribution="{params.db}/database{params.read_length}mers.kmer_distrib"
          if [ ! -s "$kmer_distribution" ]; then
            echo "Missing Bracken k-mer distribution file: $kmer_distribution" > {log}
            exit 1
          fi
          est_abundance.py \
            -i {input.report} \
            -k "$kmer_distribution" \
            -o {output.abundance} \
            -l {params.level} \
            -t {params.threshold} > {log} 2>&1
          cp {output.abundance} {output.report}
        else
          echo "Neither bracken nor est_abundance.py is available in the active environment." > {log}
          exit 127
        fi
        """


rule bracken_genus:
    input:
        report="results/taxonomy/kraken2/{sample}.kraken2.report.txt",
    output:
        abundance="results/taxonomy/bracken/{sample}.genus.bracken.tsv",
        report="results/taxonomy/bracken/{sample}.genus.kraken2.report.txt",
    log:
        "results/logs/bracken/{sample}.genus.log",
    threads:
        config["resources"].get("bracken_threads", 1)
    params:
        db=config["bracken"].get("database") or config["kraken2"].get("database"),
        read_length=config["bracken"].get("read_length", 150),
        threshold=config["bracken"].get("threshold", 10),
        level=config["bracken"].get("genus_level", "G"),
        extra=config["bracken"].get("extra_args", ""),
    conda:
        "../envs/kraken_bracken.yaml"
    shell:
        """
        mkdir -p results/taxonomy/bracken results/logs/bracken
        if [ -z "{params.db}" ] || [ "{params.db}" = "None" ] || [ "{params.db}" = "null" ]; then
          echo "Bracken database is not configured. Set bracken.database in config/config.yaml." > {log}
          exit 1
        fi
        if command -v bracken >/dev/null 2>&1; then
          bracken \
            -d {params.db} \
            -i {input.report} \
            -o {output.abundance} \
            -w {output.report} \
            -r {params.read_length} \
            -l {params.level} \
            -t {params.threshold} \
            {params.extra} > {log} 2>&1
        elif command -v est_abundance.py >/dev/null 2>&1; then
          kmer_distribution="{params.db}/database{params.read_length}mers.kmer_distrib"
          if [ ! -s "$kmer_distribution" ]; then
            echo "Missing Bracken k-mer distribution file: $kmer_distribution" > {log}
            exit 1
          fi
          est_abundance.py \
            -i {input.report} \
            -k "$kmer_distribution" \
            -o {output.abundance} \
            -l {params.level} \
            -t {params.threshold} > {log} 2>&1
          cp {output.abundance} {output.report}
        else
          echo "Neither bracken nor est_abundance.py is available in the active environment." > {log}
          exit 127
        fi
        """


rule collect_read_qc_taxonomy:
    input:
        samplesheet=config["paths"]["samplesheet"],
        fastp_json=expand("results/qc/fastp/{sample}.fastp.json", sample=SAMPLES),
        kraken_reports=expand("results/taxonomy/kraken2/{sample}.kraken2.report.txt", sample=SAMPLES),
        bracken_species=expand("results/taxonomy/bracken/{sample}.species.bracken.tsv", sample=SAMPLES),
        bracken_genus=expand("results/taxonomy/bracken/{sample}.genus.bracken.tsv", sample=SAMPLES),
    output:
        summary="results/tables/read_qc_taxonomy_summary.tsv",
        species_long="results/tables/bracken_species_abundance_long.tsv",
    log:
        "results/logs/collect_read_qc_taxonomy.log",
    params:
        thresholds="config/qc_thresholds.yaml",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/collect_kraken_bracken.py \
          --samplesheet {input.samplesheet} \
          --fastp-dir results/qc/fastp \
          --kraken-dir results/taxonomy/kraken2 \
          --bracken-dir results/taxonomy/bracken \
          --thresholds {params.thresholds} \
          --summary-output {output.summary} \
          --species-long-output {output.species_long} > {log} 2>&1
        """


rule classify_read_level_taxonomy:
    input:
        summary="results/tables/read_qc_taxonomy_summary.tsv",
    output:
        judgement="results/tables/read_level_taxonomy_initial_judgement.tsv",
    log:
        "results/logs/classify_read_level_taxonomy.log",
    params:
        thresholds="config/qc_thresholds.yaml",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/classify_read_level_taxonomy.py \
          --summary {input.summary} \
          --thresholds {params.thresholds} \
          --output {output.judgement} > {log} 2>&1
        """
