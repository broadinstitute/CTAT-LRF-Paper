#!/usr/bin/env python3

import sys, os, re
import gzip
import csv


gene_to_chr_coords = dict()
# parse chromosome coords
with gzip.open(
    "../../../../../LR-FusionBenchmarking/resources/genes.coords.gz", "rt"
) as fh:
    header = next(fh)
    for line in fh:
        line = line.rstrip()
        gene_id, chrom, lend, rend, file, primary_target = line.split("\t")
        if primary_target == "1" or gene_id not in gene_to_chr_coords:
            chrom_coords = f"{chrom}:{lend}-{rend}"
            gene_to_chr_coords[gene_id] = chrom_coords


header = ""
all_fusions = list()
with open("../Illumina_supported_fusions.tsv", "rt") as fh:
    csv_reader = csv.DictReader(fh, delimiter="\t")
    column_headers = csv_reader.fieldnames
    column_headers.extend(
        ["left_gene_coords", "right_gene_coords", "priority_fusion_name"]
    )

    writer = csv.DictWriter(
        sys.stdout, fieldnames=column_headers, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()

    for row in csv_reader:

        lex_sorted_fusion_name = row["lex_ordered_fusion_name"]

        sample_name, fusion_name = lex_sorted_fusion_name.split("|")

        left_gene, right_gene = fusion_name.split("--")

        left_gene_coords = (
            gene_to_chr_coords[left_gene] if left_gene in gene_to_chr_coords else "NA"
        )
        right_gene_coords = (
            gene_to_chr_coords[right_gene] if right_gene in gene_to_chr_coords else "NA"
        )

        row["left_gene_coords"] = left_gene_coords
        row["right_gene_coords"] = right_gene_coords

        starF_fusion_name = row["FusionName.starF"]
        arriba_fusion_name = row["FusionName.arriba"]

        priority_fusion_name = (
            starF_fusion_name if starF_fusion_name != "NA" else arriba_fusion_name
        )

        row["priority_fusion_name"] = priority_fusion_name

        writer.writerow(row)


sys.exit(0)
