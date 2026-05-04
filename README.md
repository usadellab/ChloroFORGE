<div align="left">
<img width="400" height="400" alt="ChloroFORGE_logo" src="https://github.com/user-attachments/assets/67d5e396-2b4a-42b6-8695-95ea13f169dd" />

  **Reassemble chloroplast genomes from precomputed assemblies**
</div>

ChloroFORGE is a bioinformatics pipeline that identifies chloroplast-derived contigs from existing genome assemblies (e.g., generated with [hifiasm](https://github.com/chhylp123/hifiasm) or [verkko](https://github.com/marbl/verkko) using long reads from Oxford Nanopore or PacBio) and reconstructs the complete chloroplast genome using Flye.

---

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
  - [Required Arguments](#required-arguments)
  - [Optional Arguments](#optional-arguments)
- [Examples](#examples)
- [Output Structure](#output-structure)
- [Notes & Tips](#notes--tips)
- [Citation](#citation)

---


## Dependencies

ChloroFORGE relies on the following tools (installed automatically via `setup.sh`):

| Tool | Version | Link |
|------|---------|------|
| `minimap2` | **2.30** | https://github.com/lh3/minimap2 |
| `Flye` | **2.9.6**| https://github.com/mikolmogorov/Flye |
| `seqkit` | **2.13.0** | https://github.com/shenwei356/seqkit |
| `blastn` (BLAST+) | **2.17.0** | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| `python3` | >=3.7 | — |

---

## Installation

```bash
git clone https://github.com/usadellab/chloroFORGE.git
cd chloroFORGE
bash setup.sh
```

The setup script will verify or install all required dependencies and prepare the internal pipeline structure.

---

## Usage

```bash
./chloroFORGE.sh -g GENOME -c CHLOROPLAST -t THREADS -o OUTPUT [OPTIONS]
```

### Required Arguments

| Flag | Description |
|------|-------------|
| `-o` | Sample name / project identifier |
| `-g` | Absolute path to genome assembly (FASTA) |
| `-c` | Absolute path to chloroplast reference genome (FASTA) |
| `-t` | Number of threads |

> You can download a chloroplast reference from a closely related species via [NCBI GenBank or RefSeq](https://www.ncbi.nlm.nih.gov/). Search for your organism and download the chloroplast genome in FASTA format. A useful search strategy is to combine the latin organism name with terms such as "chloroplast" and "complete genome". For example, the following query can help identify suitable reference sequences: 
  `("replace_with_species"[Organism]) AND chloroplast[All Fields] AND "complete genome"[Title] `
Use this chloroplast genome as an input for `-c`.

### Optional Arguments

| Flag | Description | Default |
|------|-------------|---------|
| `-s` | Estimated chloroplast genome size | `150k` |
| `-l` | List of contigs representing chromosomes | none |
| `-x` | Target chloroplast contig coverage for Flye | `50` |
| `-f` | Minimum contig overlap for Flye assembly | `5000` |
| `--allow-lowcov` | Allow assembly even if coverage is below target | `false` |

---

## Examples

### Basic run

```bash
./chloroFORGE.sh \
  -g genome.fasta \
  -c chloroplast_ref.fasta \
  -t 16 \
  -o sample01
```

---

### High coverage + custom overlap

For samples with large unanchored contigs, increasing coverage (`-x`) and adjusting the minimum overlap (`-f`) can improve assembly quality.

```bash
./chloroFORGE.sh \
  -g genome.fasta \
  -c chloroplast_ref.fasta \
  -t 16 \
  -o sample01 \
  -x 80 \
  -f 7000
```

---

### Low coverage assembly

If coverage is consistently low, you can either reduce `-x` and rerun, or use `--allow-lowcov` to bypass the coverage check entirely.

> **Warning:** Results from low-coverage assemblies should be inspected carefully, as they may be incomplete or incorrect.

```bash
./chloroFORGE.sh \
  -g genome.fasta \
  -c chloroplast_ref.fasta \
  -t 16 \
  -o sample01 \
  --allow-lowcov
```

---

## Output Structure

After a successful run, the output directory will have the following structure:

```
sample01/
├── blast_hits.tsv                          # Raw BLAST hits
├── cp_contigs.txt                          # Contig IDs identified as chloroplastic
├── cp_contigs.fasta                        # Extracted chloroplast contigs
├── flye_cp_out/                            # Flye output directory
├── results_cp/
│   ├── edges.fa                            # Flye assembly graph edges
│   ├── edges_depth                         # Per-edge depth file
│   └── chloroplast_final_assembly/
│       └── sample01_chloroplast.fasta      # Final oriented chloroplast assembly
└── final_genome_sample01.fasta             # Final genome (non-cp sequences + chloroplast)
```

---


## Notes & Tips

- **Assembly validation:** The pipeline stops if the final assembly length deviates more than ±10% from the reference chloroplast length. This is a sanity check — inspect your inputs if this triggers.
- **Coverage tuning:** The default target coverage of `50x` works well in most cases. If assembly fails, try lowering `-x` before resorting to `--allow-lowcov`.
- **Overlap tuning:** The default minimum overlap (`-f 5000`) can be decreased for fragmented inputs or increased for better-resolved assemblies.
- **Reference genome:** Always use a chloroplast reference from a closely related species for best BLAST sensitivity.

---

## Citation

If you use ChloroFORGE in your research, please cite this repository:

```
ChloroFORGE – https://github.com/usadellab/ChloroFORGE
```
