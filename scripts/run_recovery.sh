#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "usage: $0 SAMPLE R1 R2 MAC_BAIT_INDEX MAC_BAIT_FASTA OUTDIR THREADS" >&2
  exit 2
fi

sample=$1
r1=$2
r2=$3
bait_index=$4
bait_fasta=$5
outdir=$6
threads=$7
script_dir=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$outdir/strict" "$outdir/meta" "$outdir/mixture"

strict_r1="$outdir/strict/${sample}.R1.fastq.gz"
strict_r2="$outdir/strict/${sample}.R2.fastq.gz"
bowtie2 --very-sensitive --no-mixed --no-discordant --no-unal -X 1000 \
  -p "$threads" -x "$bait_index" -1 "$r1" -2 "$r2" \
  | samtools view -u -f 2 -F 2304 - \
  | samtools collate -u -O - \
  | samtools fastq -n -1 "$strict_r1" -2 "$strict_r2" -0 /dev/null -s /dev/null -

spades.py --isolate --only-assembler -1 "$strict_r1" -2 "$strict_r2" \
  -t "$threads" -o "$outdir/strict/spades"
python "$script_dir/filter_assembly.py" \
  --input "$outdir/strict/spades/contigs.fasta" \
  --output "$outdir/strict/${sample}.strict.min500.fasta" \
  --summary "$outdir/strict/${sample}.assembly.tsv" \
  --sample "$sample" --route strict --min-length 500

metaspades.py --only-assembler -1 "$r1" -2 "$r2" \
  -t "$threads" -o "$outdir/meta/metaspades"
minimap2 -x asm10 --secondary=no -c -t "$threads" \
  "$bait_fasta" "$outdir/meta/metaspades/contigs.fasta" \
  > "$outdir/meta/${sample}.paf"
python "$script_dir/select_bait_contigs.py" \
  --assembly "$outdir/meta/metaspades/contigs.fasta" \
  --paf "$outdir/meta/${sample}.paf" \
  --output "$outdir/meta/${sample}.meta_mac_bin.fasta" \
  --manifest "$outdir/meta/${sample}.meta_mac_bin.tsv" \
  --min-identity 0.90 --min-aligned-fraction 0.25 --min-length 500

python "$script_dir/cap_paired_fastq.py" \
  --r1 "$r1" --r2 "$r2" --max-pairs 2000000 \
  --output-r1 "$outdir/mixture/${sample}.R1.first_2m.fastq.gz" \
  --output-r2 "$outdir/mixture/${sample}.R2.first_2m.fastq.gz"

measure_burden() {
  route=$1
  assembly=$2
  prefix="$outdir/mixture/${sample}.${route}"
  samtools faidx "$assembly"
  bowtie2-build "$assembly" "$prefix"
  bowtie2 --very-sensitive --no-unal -p "$threads" -x "$prefix" \
    -1 "$outdir/mixture/${sample}.R1.first_2m.fastq.gz" \
    -2 "$outdir/mixture/${sample}.R2.first_2m.fastq.gz" \
    | samtools sort -o "$prefix.bam" -
  samtools index "$prefix.bam"
  samtools mpileup -aa -q 20 -Q 20 -d 100000 -f "$assembly" "$prefix.bam" \
    | python "$script_dir/residual_mixture.py" \
        --sample "$sample" --route "$route" --output "$prefix.tsv"
}

measure_burden strict "$outdir/strict/${sample}.strict.min500.fasta"
measure_burden meta "$outdir/meta/${sample}.meta_mac_bin.fasta"

fastANI -q "$outdir/meta/${sample}.meta_mac_bin.fasta" \
  -r "$outdir/strict/${sample}.strict.min500.fasta" \
  -o "$outdir/${sample}.meta_to_strict.tsv"
