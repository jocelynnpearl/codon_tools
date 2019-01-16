# Codon Tools
Amino acid reverse translation and DNA optimization tool based on speciies-specific codon-use distributions.
Species-specifc data can be found on the [Codon Usage Database](http://www.kazusa.or.jp) using the [NCBI Taxonomy database](http://www.ncbi.nlm.nih.gov/taxonomy) id (e.g. 413997) or the organism's Latin name (e.g. _Escherichia coli_ B). Mapping species names to Taxonomy IDs can be done [here](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi).

## Use

```sh
$python codon_tools.py --input INPUT_LIST.fasta
```

To get started, create a conda environment from the `environment.yml` file.

```sh
conda env create -f environment.yml
```

example contents of `INPUT_LIST.fasta`:

```
>SEQ_1
ACDEFGHIKLMNPQRSTVWY
>SEQ_2
ACDEFGHIKLMNPQRSTVWY
```

## Features
1. Reverse translates input amino acid sequence to DNA.
2. Calculates the host's per-AA codon usage profile -- codons used less than a specified threshold (defaults to 10%) are dropped.
3. Compares the reverse-translated DNA sequence to the host profile, determines which codons are overused/underused.
4. Stochastically mutates codons according to host profile.
5. Processes DNA to remove unwanted features:
    * high GC content within a sliding window and across the entire sequence
    * unwanted restriction sites
    * alternate start positions (GA-rich regions 18 bp upstream of ATG/GTG/TTG)
    * 3-consecutive identical codons and 9-mer repeat chunks
    * areas with more than 4 (variable) consecutive identical bps ("local homopolymers")

The process is repeated from step 3 for a specified number of cycles (defaults to 1000) OR until the per-AA codon profile of current DNA and host profile matches (within tolerance).

## To do
- [ ] remove RNA structure from sequence
- [ ] remove predicted splice sites
- [ ] store "best" sequence by profile similarity to host and return that to the user
