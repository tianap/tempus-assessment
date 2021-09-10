# Tempus Coding Challenge  - Variant Annotation

Link to repo: https://github.com/tianap/tempus-assessment

## Description
This Python script takes in a VCF file and annotates its variants with the following specifications:

Each variant must be annotated with the following pieces of information:
1. Type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent,
intergenic, etc.). If there are multiple effects, annotate with the most deleterious
possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from ExAC API (API documentation is available here:
http://exac.hms.harvard.edu/).
6. Any additional annotations that you feel might be relevant.

## Usage
```py annotateVariants.py -i '.\Challenge_data_(1).vcf' -o .\results.tsv```


