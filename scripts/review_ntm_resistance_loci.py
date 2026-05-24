#!/usr/bin/env python3
"""Create manual-review NTM resistance locus tables from local assemblies.

This module is intentionally conservative. It reports sequence evidence for
known NTM/MAC resistance-associated loci, but it does not make clinical
resistance calls without AST and species-specific review.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


RRNA_DBS = {
    "Mycobacterium_avium": "Mycobacterium avium",
    "Mycobacterium_intracellulare": "Mycobacterium intracellulare",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--assembly-dir", required=True, type=Path)
    parser.add_argument("--ntm-profiler-db-dir", required=True, type=Path)
    parser.add_argument("--priority-gff-dir", required=True, type=Path)
    parser.add_argument("--amrfinder-db-dir", required=True, type=Path)
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--tables-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--blast-bin-dir", type=Path)
    parser.add_argument("--min-blastn-identity", type=float, default=85.0)
    parser.add_argument("--min-blastn-qcov", type=float, default=70.0)
    parser.add_argument("--min-erm-identity", type=float, default=80.0)
    parser.add_argument("--min-erm-qcov", type=float, default=80.0)
    return parser.parse_args()


def run_cmd(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def tool_path(tool: str, blast_bin_dir: Path | None) -> str:
    if blast_bin_dir:
        candidate = blast_bin_dir / tool
        if candidate.exists():
            return str(candidate)
    return tool


def parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in raw.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = unquote(value)
    return attrs


def sanitize_id(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:+-]+", "_", str(text)).strip("_")


def get_subseq(records: dict[str, SeqRecord], seqid: str, start_1based: int, end_1based: int, strand: str) -> Seq:
    seq = records[seqid].seq[start_1based - 1 : end_1based]
    if strand == "-":
        seq = seq.reverse_complement()
    return seq


def load_fasta_records(path: Path) -> dict[str, SeqRecord]:
    return SeqIO.to_dict(SeqIO.parse(path, "fasta"))


def extract_ntm_profiler_rrna_refs(
    db_dir: Path, reference_dir: Path
) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    records: list[SeqRecord] = []
    meta_rows: list[dict[str, Any]] = []
    hotspot_rows: list[dict[str, Any]] = []

    for db_name, species_name in RRNA_DBS.items():
        species_dir = db_dir / db_name
        genome = load_fasta_records(species_dir / "genome.fasta")
        genes = pd.read_csv(
            species_dir / "genes.bed",
            sep="\t",
            names=["chrom", "start0", "end", "gene_id", "locus", "drug"],
        )
        gff_features: dict[str, dict[str, Any]] = {}
        for line in (species_dir / "genome.gff").open():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "rRNA":
                continue
            attrs = parse_gff_attributes(fields[8])
            gene_id = attrs.get("locus_tag") or attrs.get("Parent", "").replace("gene-", "")
            product = attrs.get("product", "")
            if product not in {"16S ribosomal RNA", "23S ribosomal RNA"}:
                continue
            gff_features[gene_id] = {
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "product": product,
            }
        mutations = json.loads((species_dir / "mutations.json").read_text())
        gene_to_locus = dict(zip(genes["gene_id"], genes["locus"]))
        gene_to_drug = dict(zip(genes["gene_id"], genes["drug"]))
        for row in genes.itertuples(index=False):
            if row.locus not in {"rrl", "rrs"}:
                continue
            feature = gff_features.get(row.gene_id)
            if not feature:
                continue
            seq = get_subseq(genome, feature["chrom"], feature["start"], feature["end"], feature["strand"])
            ref_id = f"ntmdb__{db_name}__{row.locus}__{row.gene_id}"
            records.append(
                SeqRecord(
                    seq,
                    id=ref_id,
                    description=f"{species_name} {row.locus} {row.gene_id} NTM-Profiler",
                )
            )
            meta_rows.append(
                {
                    "reference_id": ref_id,
                    "reference_source": "NTM-Profiler",
                    "db_name": db_name,
                    "species": species_name,
                    "locus": row.locus,
                    "gene_id": row.gene_id,
                    "drug_context": gene_to_drug.get(row.gene_id, ""),
                    "chrom": feature["chrom"],
                    "start": feature["start"],
                    "end": feature["end"],
                    "strand": feature["strand"],
                    "product": feature["product"],
                    "length": len(seq),
                }
            )

        for gene_id, gene_mutations in mutations.items():
            locus = gene_to_locus.get(gene_id, "unknown")
            for mutation, payload in gene_mutations.items():
                match = re.match(r"n\.(\d+)([ACGTN])>([ACGTN])", mutation)
                if not match:
                    continue
                ref_pos, ref_base, alt_base = match.groups()
                annotations = payload.get("annotations", [{}]) or [{}]
                for annot in annotations:
                    hotspot_rows.append(
                        {
                            "db_name": db_name,
                            "species": species_name,
                            "gene_id": gene_id,
                            "locus": locus,
                            "ntm_profiler_mutation": mutation,
                            "reference_position_1based": int(ref_pos),
                            "reference_base": ref_base,
                            "alternate_base": alt_base,
                            "drug": annot.get("drug", ""),
                            "annotation_type": annot.get("type", ""),
                            "e_coli_nomenclature": annot.get("e.coli-nomenclature", ""),
                            "literature": annot.get("literature", ""),
                            "comment": annot.get("comment", ""),
                        }
                    )

    reference_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = reference_dir / "ntm_profiler_rrl_rrs_reference_loci.fasta"
    SeqIO.write(records, fasta_path, "fasta")
    meta = pd.DataFrame(meta_rows)
    hotspots = pd.DataFrame(hotspot_rows)
    meta.to_csv(reference_dir / "ntm_profiler_rrl_rrs_reference_metadata.tsv", sep="\t", index=False)
    hotspots.to_csv(reference_dir / "ntm_profiler_rrl_rrs_hotspot_definitions.tsv", sep="\t", index=False)
    return fasta_path, meta, hotspots


def extract_gyr_refs(priority_gff_dir: Path, reference_dir: Path) -> tuple[Path, pd.DataFrame]:
    data_dir = priority_gff_dir / "ncbi_dataset" / "data"
    records: list[SeqRecord] = []
    rows: list[dict[str, Any]] = []
    for gff_path in sorted(data_dir.glob("GCF_*/genomic.gff")):
        accession = gff_path.parent.name
        fasta_paths = sorted(gff_path.parent.glob("*_genomic.fna"))
        if not fasta_paths:
            continue
        genome = load_fasta_records(fasta_paths[0])
        copy_counts: dict[str, int] = defaultdict(int)
        for line in gff_path.open():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = parse_gff_attributes(fields[8])
            gene = attrs.get("gene", "")
            if gene not in {"gyrA", "gyrB"}:
                continue
            if attrs.get("pseudo", "").lower() == "true":
                continue
            seqid, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            copy_counts[gene] += 1
            seq = get_subseq(genome, seqid, start, end, strand)
            ref_id = f"public__{accession}__{gene}__copy{copy_counts[gene]}"
            product = attrs.get("product", "")
            protein_id = attrs.get("protein_id", attrs.get("Name", ""))
            records.append(
                SeqRecord(
                    seq,
                    id=ref_id,
                    description=f"{accession} {gene} {protein_id} {product}",
                )
            )
            rows.append(
                {
                    "reference_id": ref_id,
                    "reference_source": "priority14_public_best_hit_gff",
                    "accession": accession,
                    "locus": gene,
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": len(seq),
                    "protein_id": protein_id,
                    "product": product,
                }
            )
    fasta_path = reference_dir / "public_best_hit_gyrA_gyrB_reference_loci.fasta"
    SeqIO.write(records, fasta_path, "fasta")
    meta = pd.DataFrame(rows)
    meta.to_csv(reference_dir / "public_best_hit_gyrA_gyrB_reference_metadata.tsv", sep="\t", index=False)
    return fasta_path, meta


def extract_erm_refs(amrfinder_db_dir: Path, reference_dir: Path) -> tuple[Path, pd.DataFrame]:
    amrprot = amrfinder_db_dir / "AMRProt.fa"
    records: list[SeqRecord] = []
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for rec in SeqIO.parse(amrprot, "fasta"):
        parts = rec.id.split("|")
        gene_symbol = parts[3] if len(parts) > 3 else rec.id
        element_symbol = parts[4] if len(parts) > 4 else gene_symbol
        drug_class = parts[8] if len(parts) > 8 else ""
        description = parts[9] if len(parts) > 9 else rec.description
        if not gene_symbol.lower().startswith("erm"):
            continue
        # Keep one representative per exact symbol and sequence to limit noisy tblastn output.
        key = (gene_symbol, str(rec.seq))
        if key in seen:
            continue
        seen.add(key)
        ref_id = f"amrfinderplus__{sanitize_id(gene_symbol)}__{rec.id.split('|')[0]}"
        records.append(SeqRecord(rec.seq, id=ref_id, description=rec.description))
        rows.append(
            {
                "reference_id": ref_id,
                "reference_source": "AMRFinderPlus_AMRProt",
                "gene_symbol": gene_symbol,
                "element_symbol": element_symbol,
                "drug_class": drug_class,
                "description": description,
                "length_aa": len(rec.seq),
            }
        )
    fasta_path = reference_dir / "amrfinderplus_erm_protein_refs.faa"
    SeqIO.write(records, fasta_path, "fasta")
    meta = pd.DataFrame(rows)
    meta.to_csv(reference_dir / "amrfinderplus_erm_protein_reference_metadata.tsv", sep="\t", index=False)
    return fasta_path, meta


def combine_fastas(paths: list[Path], out: Path) -> None:
    records = []
    for path in paths:
        records.extend(list(SeqIO.parse(path, "fasta")))
    SeqIO.write(records, out, "fasta")


def write_combined_assemblies(sample_table: pd.DataFrame, assembly_dir: Path, out: Path) -> None:
    records: list[SeqRecord] = []
    for row in sample_table.itertuples(index=False):
        sample_id = str(row.sample_id)
        assembly = assembly_dir / sample_id / f"{sample_id}.assembly.fasta"
        if not assembly.exists():
            raise FileNotFoundError(assembly)
        for rec in SeqIO.parse(assembly, "fasta"):
            rec_id = f"{sample_id}|{sanitize_id(rec.id)}"
            records.append(SeqRecord(rec.seq, id=rec_id, description=""))
    SeqIO.write(records, out, "fasta")


def run_blastn(
    ref_fasta: Path,
    assembly_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    threads: int,
    blast_bin_dir: Path | None,
) -> None:
    run_cmd([tool_path("makeblastdb", blast_bin_dir), "-in", str(assembly_fasta), "-dbtype", "nucl", "-out", str(db_prefix)])
    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
        "evalue bitscore qlen slen qcovhsp"
    )
    run_cmd(
        [
            tool_path("blastn", blast_bin_dir),
            "-query",
            str(ref_fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(out_tsv),
            "-outfmt",
            outfmt,
            "-evalue",
            "1e-20",
            "-num_threads",
            str(threads),
            "-max_target_seqs",
            "100000",
        ]
    )


def run_tblastn(
    ref_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    threads: int,
    blast_bin_dir: Path | None,
) -> None:
    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
        "evalue bitscore qlen slen qcovhsp"
    )
    run_cmd(
        [
            tool_path("tblastn", blast_bin_dir),
            "-query",
            str(ref_fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(out_tsv),
            "-outfmt",
            outfmt,
            "-evalue",
            "1e-10",
            "-num_threads",
            str(threads),
            "-max_target_seqs",
            "100000",
        ]
    )


def read_blast(path: Path) -> pd.DataFrame:
    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "qcovhsp",
    ]
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path, sep="\t", names=cols)
    for col in ["pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovhsp"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["sample_id"] = df["sseqid"].astype(str).str.split("|").str[0]
    return df


def choose_best_hits(
    blast: pd.DataFrame,
    ref_meta: pd.DataFrame,
    samples: list[str],
    loci: list[str],
    min_identity: float,
    min_qcov: float,
) -> pd.DataFrame:
    if blast.empty:
        blast = pd.DataFrame(columns=["qseqid", "sample_id"])
    annotated = blast.merge(ref_meta.rename(columns={"reference_id": "qseqid"}), on="qseqid", how="left")
    annotated = annotated[
        annotated["pident"].fillna(0).ge(min_identity) & annotated["qcovhsp"].fillna(0).ge(min_qcov)
    ].copy()
    rows: list[dict[str, Any]] = []
    for sample_id in samples:
        for locus in loci:
            sub = annotated[(annotated["sample_id"] == sample_id) & (annotated["locus"] == locus)].copy()
            if sub.empty:
                rows.append({"sample_id": sample_id, "locus": locus, "hit_status": "missing_or_below_threshold"})
                continue
            sub = sub.sort_values(["bitscore", "qcovhsp", "pident"], ascending=[False, False, False])
            row = sub.iloc[0].to_dict()
            row["hit_status"] = "hit"
            rows.append(row)
    return pd.DataFrame(rows)


def extract_hit_sequence(combined_records: dict[str, SeqRecord], hit: pd.Series) -> Seq:
    contig = str(hit["sseqid"])
    start = int(hit["sstart"])
    end = int(hit["send"])
    lo, hi = min(start, end), max(start, end)
    seq = combined_records[contig].seq[lo - 1 : hi]
    if start > end:
        seq = seq.reverse_complement()
    return seq


def make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def alignment_position_map(ref_seq: Seq, query_seq: Seq) -> dict[int, str]:
    aligner = make_aligner()
    alignment = aligner.align(str(ref_seq).upper(), str(query_seq).upper())[0]
    ref_aln = alignment[0, :]
    query_aln = alignment[1, :]
    pos_map: dict[int, str] = {}
    ref_pos = 0
    for ref_base, query_base in zip(ref_aln, query_aln):
        if ref_base != "-":
            ref_pos += 1
            pos_map[ref_pos] = query_base
    return pos_map


def all_aligned_bases(ref_seq: Seq, query_seq: Seq) -> tuple[str, str]:
    aligner = make_aligner()
    alignment = aligner.align(str(ref_seq).upper(), str(query_seq).upper())[0]
    return alignment[0, :], alignment[1, :]


def assign_hotspot_db(species: str) -> tuple[str, str]:
    if species == "Mycobacterium avium":
        return "Mycobacterium_avium", "species_specific"
    return "Mycobacterium_intracellulare", "MAC_near_neighbor_db_caution"


def call_hotspots(
    sample_table: pd.DataFrame,
    blast: pd.DataFrame,
    rrna_meta: pd.DataFrame,
    hotspots: pd.DataFrame,
    combined_records: dict[str, SeqRecord],
    ref_records: dict[str, SeqRecord],
    min_identity: float,
    min_qcov: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    blast = blast.merge(rrna_meta.rename(columns={"reference_id": "qseqid"}), on="qseqid", how="left")
    for sample in sample_table.itertuples(index=False):
        sample_id = str(sample.sample_id)
        species = str(sample.amr_module_species)
        db_name, db_note = assign_hotspot_db(species)
        for locus in ["rrl", "rrs"]:
            ref_row = rrna_meta[(rrna_meta["db_name"] == db_name) & (rrna_meta["locus"] == locus)]
            if ref_row.empty:
                continue
            ref_id = ref_row.iloc[0]["reference_id"]
            sub = blast[
                (blast["sample_id"] == sample_id)
                & (blast["qseqid"] == ref_id)
                & blast["pident"].fillna(0).ge(min_identity)
                & blast["qcovhsp"].fillna(0).ge(min_qcov)
            ].copy()
            locus_hotspots = hotspots[(hotspots["db_name"] == db_name) & (hotspots["locus"] == locus)]
            if sub.empty:
                for hs in locus_hotspots.itertuples(index=False):
                    rows.append(
                        {
                            "sample_id": sample_id,
                            "species": species,
                            "hotspot_db": db_name,
                            "hotspot_db_note": db_note,
                            "locus": locus,
                            "drug": hs.drug,
                            "ntm_profiler_mutation": hs.ntm_profiler_mutation,
                            "e_coli_nomenclature": hs.e_coli_nomenclature,
                            "reference_position_1based": hs.reference_position_1based,
                            "reference_base": hs.reference_base,
                            "alternate_base": hs.alternate_base,
                            "sample_base": "NA",
                            "hotspot_status": "locus_missing_or_below_threshold",
                            "manual_review_note": "No adequate BLASTN hit to the assigned NTM-Profiler locus reference.",
                        }
                    )
                continue
            best = sub.sort_values(["bitscore", "qcovhsp", "pident"], ascending=[False, False, False]).iloc[0]
            sample_seq = extract_hit_sequence(combined_records, best)
            ref_seq = ref_records[ref_id].seq
            pos_map = alignment_position_map(ref_seq, sample_seq)
            for hs in locus_hotspots.itertuples(index=False):
                sample_base = pos_map.get(int(hs.reference_position_1based), "NA").upper()
                if sample_base == hs.alternate_base:
                    status = "known_resistance_associated_alt_detected"
                elif sample_base == hs.reference_base:
                    status = "wildtype_at_this_hotspot"
                elif sample_base in {"A", "C", "G", "T"}:
                    status = "other_base_at_hotspot"
                else:
                    status = "uncalled_or_gap_at_hotspot"
                rows.append(
                    {
                        "sample_id": sample_id,
                        "species": species,
                        "hotspot_db": db_name,
                        "hotspot_db_note": db_note,
                        "locus": locus,
                        "drug": hs.drug,
                        "ntm_profiler_mutation": hs.ntm_profiler_mutation,
                        "e_coli_nomenclature": hs.e_coli_nomenclature,
                        "reference_position_1based": hs.reference_position_1based,
                        "reference_base": hs.reference_base,
                        "alternate_base": hs.alternate_base,
                        "sample_base": sample_base,
                        "hotspot_status": status,
                        "pident": best["pident"],
                        "qcovhsp": best["qcovhsp"],
                        "sseqid": best["sseqid"],
                        "sstart": best["sstart"],
                        "send": best["send"],
                        "manual_review_note": (
                            "Species-specific NTM-Profiler reference used."
                            if db_note == "species_specific"
                            else "Nearest available MAC NTM-Profiler reference used; interpret species-specific effect cautiously."
                        ),
                    }
                )
    return pd.DataFrame(rows)


def call_coding_variants(
    sample_table: pd.DataFrame,
    best_hits: pd.DataFrame,
    combined_records: dict[str, SeqRecord],
    ref_records: dict[str, SeqRecord],
    min_interpretable_identity: float = 95.0,
    min_interpretable_qcov: float = 90.0,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    species_map = dict(zip(sample_table["sample_id"], sample_table["amr_module_species"]))
    for hit in best_hits.itertuples(index=False):
        hit_dict = hit._asdict()
        sample_id = str(hit_dict["sample_id"])
        locus = str(hit_dict["locus"])
        if hit_dict.get("hit_status") != "hit":
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": species_map.get(sample_id, ""),
                    "locus": locus,
                    "variant_status": "locus_missing_or_below_threshold",
                    "manual_review_note": "No adequate BLASTN hit for this coding locus.",
                }
            )
            continue
        ref_id = str(hit_dict["qseqid"])
        if float(hit_dict.get("pident") or 0) < min_interpretable_identity or float(hit_dict.get("qcovhsp") or 0) < min_interpretable_qcov:
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": species_map.get(sample_id, ""),
                    "locus": locus,
                    "reference_id": ref_id,
                    "variant_status": "low_identity_reference_for_manual_review_only",
                    "pident": hit_dict.get("pident"),
                    "qcovhsp": hit_dict.get("qcovhsp"),
                    "sseqid": hit_dict.get("sseqid"),
                    "sstart": hit_dict.get("sstart"),
                    "send": hit_dict.get("send"),
                    "manual_review_note": (
                        f"Best reference identity/coverage below interpretable threshold "
                        f"({min_interpretable_identity}% identity, {min_interpretable_qcov}% coverage); "
                        "do not count background divergence as resistance-associated mutation."
                    ),
                }
            )
            continue
        ref_seq = ref_records[ref_id].seq
        sample_seq = extract_hit_sequence(combined_records, pd.Series(hit_dict))
        ref_aln, query_aln = all_aligned_bases(ref_seq, sample_seq)
        pos_to_query: dict[int, str] = {}
        ref_pos = 0
        for ref_base, query_base in zip(ref_aln, query_aln):
            if ref_base != "-":
                ref_pos += 1
                pos_to_query[ref_pos] = query_base.upper()
        changes = 0
        codon_count = len(ref_seq) // 3
        for codon_idx in range(codon_count):
            ref_start = codon_idx * 3 + 1
            ref_codon = str(ref_seq[ref_start - 1 : ref_start + 2]).upper()
            query_codon = "".join(pos_to_query.get(pos, "-") for pos in range(ref_start, ref_start + 3)).upper()
            if query_codon == ref_codon:
                continue
            if "-" in query_codon or "N" in query_codon or len(query_codon) != 3:
                ref_aa = str(Seq(ref_codon).translate())
                query_aa = "X"
                effect = "uncalled_or_indel_codon"
            else:
                ref_aa = str(Seq(ref_codon).translate())
                query_aa = str(Seq(query_codon).translate())
                effect = "synonymous" if ref_aa == query_aa else "nonsynonymous"
            if effect == "synonymous":
                continue
            changes += 1
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": species_map.get(sample_id, ""),
                    "locus": locus,
                    "reference_id": ref_id,
                    "reference_codon_position": codon_idx + 1,
                    "reference_nt_position_start_1based": ref_start,
                    "reference_codon": ref_codon,
                    "sample_codon": query_codon,
                    "reference_aa": ref_aa,
                    "sample_aa": query_aa,
                    "aa_change": f"{ref_aa}{codon_idx + 1}{query_aa}",
                    "variant_status": effect,
                    "pident": hit_dict.get("pident"),
                    "qcovhsp": hit_dict.get("qcovhsp"),
                    "sseqid": hit_dict.get("sseqid"),
                    "sstart": hit_dict.get("sstart"),
                    "send": hit_dict.get("send"),
                    "manual_review_note": "Amino-acid change relative to nearest public/reference coding locus; not a clinical resistance call.",
                }
            )
        if changes == 0:
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": species_map.get(sample_id, ""),
                    "locus": locus,
                    "reference_id": ref_id,
                    "variant_status": "no_nonsynonymous_changes_detected",
                    "pident": hit_dict.get("pident"),
                    "qcovhsp": hit_dict.get("qcovhsp"),
                    "sseqid": hit_dict.get("sseqid"),
                    "sstart": hit_dict.get("sstart"),
                    "send": hit_dict.get("send"),
                    "manual_review_note": "No nonsynonymous coding changes relative to the selected nearest reference locus.",
                }
            )
    return pd.DataFrame(rows)


def summarize_erm(
    sample_table: pd.DataFrame,
    tblastn: pd.DataFrame,
    erm_meta: pd.DataFrame,
    min_identity: float,
    min_qcov: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    annotated = tblastn.merge(erm_meta.rename(columns={"reference_id": "qseqid"}), on="qseqid", how="left")
    for sample in sample_table.itertuples(index=False):
        sample_id = str(sample.sample_id)
        sub = annotated[
            (annotated["sample_id"] == sample_id)
            & annotated["pident"].fillna(0).ge(min_identity)
            & annotated["qcovhsp"].fillna(0).ge(min_qcov)
        ].copy()
        if sub.empty:
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": sample.amr_module_species,
                    "erm_screen_status": "no_confident_erm_hit",
                    "manual_review_note": f"No tblastn hit meeting identity>={min_identity}% and query_coverage>={min_qcov}%.",
                }
            )
            continue
        sub = sub.sort_values(["bitscore", "qcovhsp", "pident"], ascending=[False, False, False])
        for _, row in sub.head(5).iterrows():
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": sample.amr_module_species,
                    "erm_screen_status": "confident_erm_like_hit",
                    "reference_id": row["qseqid"],
                    "gene_symbol": row.get("gene_symbol", ""),
                    "element_symbol": row.get("element_symbol", ""),
                    "description": row.get("description", ""),
                    "pident": row["pident"],
                    "qcovhsp": row["qcovhsp"],
                    "bitscore": row["bitscore"],
                    "sseqid": row["sseqid"],
                    "sstart": row["sstart"],
                    "send": row["send"],
                    "manual_review_note": "tblastn screen against AMRFinderPlus Erm protein references; review for specificity before interpretation.",
                }
            )
    return pd.DataFrame(rows)


def build_manual_summary(
    sample_table: pd.DataFrame,
    hotspot_table: pd.DataFrame,
    coding_table: pd.DataFrame,
    erm_table: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for sample in sample_table.itertuples(index=False):
        sample_id = str(sample.sample_id)
        hs = hotspot_table[hotspot_table["sample_id"] == sample_id]
        known = hs[hs["hotspot_status"].eq("known_resistance_associated_alt_detected")]
        uncalled = hs[hs["hotspot_status"].str.contains("uncalled|missing", na=False)]
        coding = coding_table[coding_table["sample_id"] == sample_id]
        gyr_changes = coding[coding["variant_status"].eq("nonsynonymous")]
        low_identity_gyr = coding[coding["variant_status"].eq("low_identity_reference_for_manual_review_only")]
        erm = erm_table[erm_table["sample_id"] == sample_id]
        erm_conf = erm["erm_screen_status"].eq("confident_erm_like_hit").any()
        review_flags = []
        if not known.empty:
            review_flags.append("known_rrl_rrs_resistance_hotspot_alt")
        if not uncalled.empty:
            review_flags.append("uncalled_rrl_rrs_hotspot")
        if not gyr_changes.empty:
            review_flags.append("gyrA_gyrB_nonsynonymous_change")
        if not low_identity_gyr.empty:
            review_flags.append("low_identity_gyr_reference_review_only")
        if erm_conf:
            review_flags.append("erm_like_hit")
        if not review_flags:
            review_flags.append("no_known_hotspot_alt_detected")
        rows.append(
            {
                "sample_id": sample_id,
                "species": sample.amr_module_species,
                "known_rrl_rrs_hotspot_alt_count": len(known),
                "uncalled_or_missing_hotspot_count": len(uncalled),
                "gyrA_gyrB_nonsynonymous_change_count": len(gyr_changes),
                "low_identity_gyr_reference_count": len(low_identity_gyr),
                "confident_erm_like_hit": bool(erm_conf),
                "manual_review_flags": ";".join(review_flags),
                "interpretation_boundary": "Genomic review table only; do not infer clinical resistance without AST/species-specific validation.",
            }
        )
    return pd.DataFrame(rows)


def build_mutation_review_table(
    sample_table: pd.DataFrame,
    hotspot_table: pd.DataFrame,
    coding_table: pd.DataFrame,
    erm_table: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    species_map = dict(zip(sample_table["sample_id"], sample_table["amr_module_species"]))
    boundary = "Manual genomic review only; not a clinical resistance call without AST/MIC and species-specific validation."

    for sample_id in sample_table["sample_id"].astype(str):
        species = species_map.get(sample_id, "")
        sample_hotspots = hotspot_table[hotspot_table["sample_id"] == sample_id]
        for locus in ["rrl", "rrs"]:
            sub = sample_hotspots[sample_hotspots["locus"] == locus].copy()
            drug_context = "macrolides" if locus == "rrl" else "aminoglycosides"
            notable = sub[~sub["hotspot_status"].eq("wildtype_at_this_hotspot")].copy()
            if notable.empty:
                ref_system = "NA"
                if not sub.empty:
                    ref_system = f"{sub.iloc[0]['hotspot_db']} ({sub.iloc[0]['hotspot_db_note']})"
                rows.append(
                    {
                        "sample_id": sample_id,
                        "species": species,
                        "locus": locus,
                        "drug_or_class_context": drug_context,
                        "evidence_layer": "NTM-Profiler curated rRNA hotspot",
                        "review_status": "all_defined_hotspots_wildtype",
                        "mutation_or_change": "all_defined_NTM-Profiler_hotspots_wildtype",
                        "reference_system": ref_system,
                        "reference_id": "NTM-Profiler mutations.json",
                        "pident": sub["pident"].max() if "pident" in sub else pd.NA,
                        "qcovhsp": sub["qcovhsp"].max() if "qcovhsp" in sub else pd.NA,
                        "contig_or_subject": sub.iloc[0]["sseqid"] if not sub.empty and "sseqid" in sub else "NA",
                        "coordinate_hint": "defined hotspot set",
                        "manual_review_note": "No alternate base was detected at the curated hotspot definitions available in the local NTM-Profiler database.",
                        "interpretation_boundary": boundary,
                    }
                )
            else:
                for _, row in notable.iterrows():
                    rows.append(
                        {
                            "sample_id": sample_id,
                            "species": species,
                            "locus": locus,
                            "drug_or_class_context": row.get("drug", drug_context),
                            "evidence_layer": "NTM-Profiler curated rRNA hotspot",
                            "review_status": row.get("hotspot_status", ""),
                            "mutation_or_change": row.get("e_coli_nomenclature") or row.get("ntm_profiler_mutation", ""),
                            "reference_system": f"{row.get('hotspot_db', '')} ({row.get('hotspot_db_note', '')})",
                            "reference_id": "NTM-Profiler mutations.json",
                            "pident": row.get("pident", pd.NA),
                            "qcovhsp": row.get("qcovhsp", pd.NA),
                            "contig_or_subject": row.get("sseqid", "NA"),
                            "coordinate_hint": row.get("reference_position_1based", ""),
                            "manual_review_note": row.get("manual_review_note", ""),
                            "interpretation_boundary": boundary,
                        }
                    )

        sample_coding = coding_table[coding_table["sample_id"] == sample_id]
        for locus in ["gyrA", "gyrB"]:
            sub = sample_coding[sample_coding["locus"] == locus].copy()
            if sub.empty:
                rows.append(
                    {
                        "sample_id": sample_id,
                        "species": species,
                        "locus": locus,
                        "drug_or_class_context": "fluoroquinolone_candidate",
                        "evidence_layer": "coding-locus reference alignment",
                        "review_status": "locus_not_evaluated",
                        "mutation_or_change": "NA",
                        "reference_system": "priority public/reference coding loci",
                        "reference_id": "NA",
                        "manual_review_note": "No coding-locus row was generated for this sample/locus.",
                        "interpretation_boundary": boundary,
                    }
                )
                continue
            notable = sub[~sub["variant_status"].eq("no_nonsynonymous_changes_detected")].copy()
            rows_to_write = notable if not notable.empty else sub.head(1)
            for _, row in rows_to_write.iterrows():
                rows.append(
                    {
                        "sample_id": sample_id,
                        "species": species,
                        "locus": locus,
                        "drug_or_class_context": "fluoroquinolone_candidate",
                        "evidence_layer": "coding-locus reference alignment",
                        "review_status": row.get("variant_status", ""),
                        "mutation_or_change": row.get("aa_change") if pd.notna(row.get("aa_change", pd.NA)) else row.get("variant_status", ""),
                        "reference_system": "nearest priority public/reference coding locus",
                        "reference_id": row.get("reference_id", "NA"),
                        "pident": row.get("pident", pd.NA),
                        "qcovhsp": row.get("qcovhsp", pd.NA),
                        "contig_or_subject": row.get("sseqid", "NA"),
                        "coordinate_hint": row.get("reference_codon_position", ""),
                        "manual_review_note": row.get("manual_review_note", ""),
                        "interpretation_boundary": boundary,
                    }
                )

        sample_erm = erm_table[erm_table["sample_id"] == sample_id]
        if sample_erm.empty:
            rows.append(
                {
                    "sample_id": sample_id,
                    "species": species,
                    "locus": "erm",
                    "drug_or_class_context": "macrolide_inducible_candidate",
                    "evidence_layer": "AMRFinderPlus Erm protein tblastn screen",
                    "review_status": "not_evaluated",
                    "mutation_or_change": "NA",
                    "reference_system": "AMRFinderPlus AMRProt",
                    "reference_id": "NA",
                    "manual_review_note": "No erm screen row was generated.",
                    "interpretation_boundary": boundary,
                }
            )
        else:
            for _, row in sample_erm.iterrows():
                rows.append(
                    {
                        "sample_id": sample_id,
                        "species": species,
                        "locus": "erm",
                        "drug_or_class_context": "macrolide_inducible_candidate",
                        "evidence_layer": "AMRFinderPlus Erm protein tblastn screen",
                        "review_status": row.get("erm_screen_status", ""),
                        "mutation_or_change": row.get("gene_symbol", "no_confident_erm_hit"),
                        "reference_system": "AMRFinderPlus AMRProt",
                        "reference_id": row.get("reference_id", "NA"),
                        "pident": row.get("pident", pd.NA),
                        "qcovhsp": row.get("qcovhsp", pd.NA),
                        "contig_or_subject": row.get("sseqid", "NA"),
                        "coordinate_hint": f"{row.get('sstart', '')}-{row.get('send', '')}",
                        "manual_review_note": row.get("manual_review_note", ""),
                        "interpretation_boundary": boundary,
                    }
                )

    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.reference_dir.mkdir(parents=True, exist_ok=True)
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.tables_dir.mkdir(parents=True, exist_ok=True)

    sample_table = pd.read_csv(args.sample_table, sep="\t")
    sample_ids = sample_table["sample_id"].astype(str).tolist()

    rrna_fasta, rrna_meta, hotspots = extract_ntm_profiler_rrna_refs(args.ntm_profiler_db_dir, args.reference_dir)
    gyr_fasta, gyr_meta = extract_gyr_refs(args.priority_gff_dir, args.reference_dir)
    erm_fasta, erm_meta = extract_erm_refs(args.amrfinder_db_dir, args.reference_dir)

    nt_refs = args.reference_dir / "ntm_resistance_nt_locus_references.fasta"
    combine_fastas([rrna_fasta, gyr_fasta], nt_refs)
    nt_meta = pd.concat([rrna_meta, gyr_meta], ignore_index=True)
    nt_meta.to_csv(args.reference_dir / "ntm_resistance_nt_locus_reference_metadata.tsv", sep="\t", index=False)

    combined_assemblies = args.workdir / "high_confidence_13_combined_assemblies.fasta"
    write_combined_assemblies(sample_table, args.assembly_dir, combined_assemblies)
    db_prefix = args.workdir / "high_confidence_13_combined_assemblies_blastdb"
    blastn_tsv = args.workdir / "ntm_resistance_loci_blastn.tsv"
    tblastn_tsv = args.workdir / "erm_tblastn.tsv"
    run_blastn(nt_refs, combined_assemblies, db_prefix, blastn_tsv, args.threads, args.blast_bin_dir)
    run_tblastn(erm_fasta, db_prefix, tblastn_tsv, args.threads, args.blast_bin_dir)

    blastn = read_blast(blastn_tsv)
    tblastn = read_blast(tblastn_tsv)
    blastn.to_csv(args.tables_dir / "ntm_resistance_loci_blastn_all_hits.tsv", sep="\t", index=False)
    tblastn.to_csv(args.tables_dir / "ntm_erm_tblastn_all_hits.tsv", sep="\t", index=False)

    best_hits = choose_best_hits(
        blastn,
        nt_meta,
        sample_ids,
        ["rrl", "rrs", "gyrA", "gyrB"],
        args.min_blastn_identity,
        args.min_blastn_qcov,
    )
    best_hits.to_csv(args.tables_dir / "ntm_resistance_locus_hit_summary.tsv", sep="\t", index=False)

    combined_records = load_fasta_records(combined_assemblies)
    rrna_records = load_fasta_records(rrna_fasta)
    gyr_records = load_fasta_records(gyr_fasta)
    nt_records = {**rrna_records, **gyr_records}

    hotspot_table = call_hotspots(
        sample_table,
        blastn,
        rrna_meta,
        hotspots,
        combined_records,
        rrna_records,
        args.min_blastn_identity,
        args.min_blastn_qcov,
    )
    hotspot_table.to_csv(args.tables_dir / "ntm_resistance_hotspot_review_table.tsv", sep="\t", index=False)

    coding_hits = best_hits[best_hits["locus"].isin(["gyrA", "gyrB"])].copy()
    coding_table = call_coding_variants(sample_table, coding_hits, combined_records, nt_records)
    coding_table.to_csv(args.tables_dir / "ntm_resistance_coding_variant_review_table.tsv", sep="\t", index=False)

    erm_table = summarize_erm(sample_table, tblastn, erm_meta, args.min_erm_identity, args.min_erm_qcov)
    erm_table.to_csv(args.tables_dir / "ntm_erm_screen_review_table.tsv", sep="\t", index=False)

    manual_summary = build_manual_summary(sample_table, hotspot_table, coding_table, erm_table)
    manual_summary.to_csv(args.tables_dir / "ntm_resistance_manual_review_summary.tsv", sep="\t", index=False)

    mutation_review = build_mutation_review_table(sample_table, hotspot_table, coding_table, erm_table)
    mutation_review.to_csv(
        args.tables_dir / "ntm_resistance_mutation_review_for_manual_curation.tsv",
        sep="\t",
        index=False,
    )

    methods = args.tables_dir / "ntm_resistance_locus_method_notes.md"
    methods.write_text(
        "\n".join(
            [
                "# NTM Resistance Locus Manual Review Notes",
                "",
                "This table set is a genomic review aid, not a clinical resistance report.",
                "",
                "- rrl/rrs known hotspot definitions were parsed from the local NTM-Profiler databases for Mycobacterium avium and Mycobacterium intracellulare.",
                "- M. avium isolates use the M. avium NTM-Profiler hotspot coordinate system.",
                "- Non-avium MAC isolates use the M. intracellulare/Chimaera NTM-Profiler hotspot coordinate system as the nearest available MAC reference and are flagged with caution.",
                "- gyrA/gyrB coding loci were extracted from already downloaded priority public-neighbor GFF3 files and compared to the nearest BLASTN reference hit.",
                "- erm was screened by tblastn against AMRFinderPlus Erm protein references with stringent identity and coverage thresholds.",
                "- ntm_resistance_mutation_review_for_manual_curation.tsv is the compact table intended for manual review; detailed hotspot and BLAST hit tables are retained for traceability.",
                "- Reported variants should be manually checked against assemblies, read support, species, and AST/MIC metadata before biological interpretation.",
                "",
            ]
        )
    )


if __name__ == "__main__":
    main()
