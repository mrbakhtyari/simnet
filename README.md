# SimNet

[![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Sequence Similarity Network Approach for Fast Detection of Horizontal Gene Transfer Events**

SimNet-HGT is a scalable computational framework for detecting horizontal gene transfer (HGT) events using sequence similarity networks (SSNs). By combining phylogenetic distance and sequence similarity with a ranking-based normalization, our method prioritizes high-similarity connections between evolutionarily distant lineages—a hallmark signature of horizontal gene transfer.

## 📖 Overview

Horizontal gene transfer is a primary driver of microbial evolution, enabling the rapid spread of crucial traits such as antibiotic resistance and virulence factors. Traditional HGT detection methods based on gene tree–species tree reconciliation suffer from prohibitive $O(nlN^3)$ time complexity, making them impractical for modern metagenomic datasets.

SimNet-HGT overcomes these limitations by:

1. **Constructing Sequence Similarity Networks** — Using MMseqs2 for rapid all-against-all sequence alignment with $O(n)$ time complexity
2. **Building Bipartite Graph Representations** — Linking genomes to connected components (gene families) for efficient modularity analysis
3. **Applying Ranking-Based Scoring** — Weighting phylogenetic distance by sequence similarity, normalized by similarity rank
4. **Statistical Validation** — Using bootstrap analysis for principled identification of high-confidence HGT candidates

The framework achieves **linear time complexity** with respect to the number of input sequences, enabling analysis of datasets comprising thousands of taxa and millions of gene/protein sequences.

## 🚀 Installation

### Prerequisites

- **Python 3.13+**
- **MMseqs2** — For sequence similarity search ([installation guide](https://github.com/soedinglab/MMseqs2#installation))

### Using uv (Recommended)

```bash
# Clone the repository
git clone https://github.com/mrbakhtyari/simnet.git
cd simnet

# Install dependencies with uv
uv sync
```

### Alternative (venv + pip)
If uv is not available, a standard virtual environment with pip also works.
```bash
python -m venv .venv
# activate the virtual environment
python -m pip install -e .
```

## 🔧 Usage
    
To run the analysis on a dataset:

1. **Prepare your dataset directory:**
   Create a root directory for your dataset (e.g., `datasets/my_dataset`).
   Inside, create a `00_input` folder containing two files:
   
   - `inputfile`: A text file containing all protein sequences.
     - **Format**:
       - Line 1: `<count> <length>`
       - Subsequent lines: `<organism_name> <protein_sequence>`
   
   - `tree.newick`: A standard Newick format phylogenetic tree containing the organisms.
     - **Note**: Organism names in the tree must match (or be substrings of) the names in `inputfile`.

2. **Run the pipeline:**

   ```bash
   python scripts/run_pipeline.py datasets/my_dataset
   ```

   **Optional Arguments:**
   - `--min_coverage`: Minimum alignment coverage
   - `--min_identity`: Minimum sequence identity

3. **Output:**
   Results will be generated in the `my_dataset` folder, organized by step:
   - `06_hgt/hgt_scores_original_names.tsv`: Final ranked HGT candidates.
   - `06_hgt/high_confidence_hgt_original_names.tsv`: High confidence subset (if applicable).

## 🧬 Execution Example: Aminoacyl

This example demonstrates how to reproduce the results for the **aminoacyl-tRNA synthetases (AARS)** dataset reported in the paper, comprising 32 organisms spanning Bacteria, Archaea, and Eukarya (Woese et al., 2000).

### Input Data
The dataset is located in `datasets/aminoacyl/00_input/`.

*   **Species Tree** (`tree.newick`):
*   **Protein Sequences** (`inputfile`):
    Contains 32 aminoacyl-tRNA synthetase sequences in Phylip format (header: `32 171`).

### Execution Command
Run the pipeline with the specific parameters used in our analysis:

```bash
python scripts/run_pipeline.py datasets/aminoacyl
```

> **Note**: This command uses the default parameters defined in the pipeline, which correspond exactly to the reported results:
> *   **Min Coverage**: 0.8
> *   **Min Identity**: 0.59
> *   **Rank Threshold**: 0.5
> *   **Bootstrap Replicates**: 10,000 (used for threshold calculation)
> *   **Target Quantile**: 95th percentile

### Expected Results

- Successfully recovered known HGT events between *Treponema pallidum*–*Pyrococcus horikoshii* and *Borrelia burgdorferi*–*Pyrococcus horikoshii*
- Identified 11 high-confidence HGT candidates exceeding the bootstrap-derived significance threshold
- Discovered novel transfer candidates among bacterial lineages

## 📝 Citation

This section will be updated upon paper acceptance.

## 👥 Authors

- **Mohammadreza Bakhtyari** — Computer Science, Université du Québec à Montréal, Canada
- **Nadia Tahiri** — Computer Science, Université de Sherbrooke, Canada
- **F. Guillaume Blanchet** — Biological Science, Université de Sherbrooke, Canada
- **Vladimir Makarenkov** — Computer Science, Université du Québec à Montréal, Canada
