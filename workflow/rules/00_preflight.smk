rule validate_samplesheet:
    input:
        samplesheet=config["paths"]["samplesheet"],
    output:
        report="results/tables/preflight_samplesheet_check.tsv",
    log:
        "results/logs/preflight.log",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/validate_samplesheet.py \
          --samplesheet {input.samplesheet} \
          --output {output.report} \
          --project-root . > {log} 2>&1
        """
