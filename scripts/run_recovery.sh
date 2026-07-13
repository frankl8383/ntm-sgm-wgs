#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 7 || $# -gt 8 ]]; then
  echo "usage: $0 SAMPLE R1 R2 MAC_BAIT_INDEX MAC_BAIT_FASTA OUTDIR THREADS [BURDEN_PAIRS]" >&2
  exit 2
fi

sample=$1
r1=$2
r2=$3
bait_index=$4
bait_fasta=$5
outdir=$6
threads=$7
burden_pairs=${8:-2000000}
script_dir=$(cd "$(dirname "$0")" && pwd)

strict_dir="$outdir/strict/$sample"
meta_dir="$outdir/meta/$sample"
mixture_dir="$outdir/mixture_metrics/$sample"
minor_input_dir="$outdir/minor_inputs/$sample"
mkdir -p "$strict_dir/logs" "$meta_dir" "$mixture_dir" "$minor_input_dir"

strict_r1="$strict_dir/${sample}.R1.fastq.gz"
strict_r2="$strict_dir/${sample}.R2.fastq.gz"
strict_fasta="$strict_dir/${sample}.strict.min500.fasta"
meta_fasta="$meta_dir/${sample}.metaspades_mac_bin.min500.fasta"
bowtie2 --very-sensitive --no-mixed --no-discordant --no-unal -X 1000 \
  -p "$threads" -x "$bait_index" -1 "$r1" -2 "$r2" \
  2> "$strict_dir/logs/${sample}.bowtie2.log" \
  | samtools view -u -f 2 -F 2304 - \
  | samtools collate -u -O - \
  | samtools fastq -n -1 "$strict_r1" -2 "$strict_r2" -0 /dev/null -s /dev/null -

spades.py --isolate --only-assembler -1 "$strict_r1" -2 "$strict_r2" \
  -t "$threads" -o "$strict_dir/spades"
python "$script_dir/filter_assembly.py" \
  --input "$strict_dir/spades/contigs.fasta" \
  --output "$strict_fasta" \
  --summary "$strict_dir/${sample}.assembly_stats.tsv" \
  --sample "$sample" --route strict --min-length 500

metaspades.py --only-assembler -1 "$r1" -2 "$r2" \
  -t "$threads" -o "$meta_dir/metaspades"
minimap2 -x asm10 --secondary=no -c -t "$threads" \
  "$bait_fasta" "$meta_dir/metaspades/contigs.fasta" \
  > "$meta_dir/${sample}.paf"
python "$script_dir/select_bait_contigs.py" \
  --assembly "$meta_dir/metaspades/contigs.fasta" \
  --paf "$meta_dir/${sample}.paf" \
  --output "$meta_dir/${sample}.metaspades_mac_bin.fasta" \
  --manifest "$meta_dir/${sample}.metaspades_mac_bin.tsv" \
  --min-identity 0.90 --min-aligned-fraction 0.25 --min-length 500
python "$script_dir/filter_assembly.py" \
  --input "$meta_dir/${sample}.metaspades_mac_bin.fasta" \
  --output "$meta_fasta" \
  --summary "$meta_dir/${sample}.metaspades_mac_bin_stats.tsv" \
  --sample "$sample" --route meta --min-length 500

python "$script_dir/cap_paired_fastq.py" \
  --r1 "$r1" --r2 "$r2" --max-pairs "$burden_pairs" \
  --output-r1 "$minor_input_dir/${sample}.R1.burden_input.fastq.gz" \
  --output-r2 "$minor_input_dir/${sample}.R2.burden_input.fastq.gz"

measure_burden() {
  route=$1
  assembly=$2
  prefix="$mixture_dir/${sample}.${route}"
  samtools faidx "$assembly"
  bowtie2-build "$assembly" "$prefix"
  bowtie2 --very-sensitive --no-unal -p "$threads" -x "$prefix" \
    -1 "$minor_input_dir/${sample}.R1.burden_input.fastq.gz" \
    -2 "$minor_input_dir/${sample}.R2.burden_input.fastq.gz" \
    | samtools sort -o "$prefix.bam" -
  samtools index "$prefix.bam"
  samtools mpileup -aa -q 20 -Q 20 -d 100000 -f "$assembly" "$prefix.bam" \
    | python "$script_dir/residual_mixture.py" \
        --sample "$sample" --route "$route" --output "$prefix.tsv"
}

measure_burden strict "$strict_fasta"
measure_burden meta "$meta_fasta"

fastANI -q "$meta_fasta" -r "$strict_fasta" \
  -o "$meta_dir/${sample}.meta_to_strict.tsv"
fastANI -q "$strict_fasta" -r "$meta_fasta" \
  -o "$meta_dir/${sample}.strict_to_meta.tsv"
