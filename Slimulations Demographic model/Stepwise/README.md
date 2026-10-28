# Stepwise Demographic Model

This folder extends the single-population pipeline at the repo root into a full
**multi-population demographic history**: a tree of populations that split off
from one another (read from a metadata table of parent/child edges, population
sizes, and generation times), each evolving short tandem repeats (STRs) under
the **stepwise mutation model (SMM)**. For every population in the tree it
extracts STR genotypes and a clean VCF, then runs the TREVA likelihood
inference to recover the mutation rate from the simulated data.

## 1. Code overview

1. **`workflow.py`**
   A `gwf` pipeline that drives the whole run for one or more mutation rates. For each mutation rate it:
   - Writes a `Makefile` that simulates every population in the demographic tree in dependency order (a child population's SLiM run takes its parent's tree sequence as input, so parents always simulate before their children).
   - Runs that `Makefile` through a single `simulate` gwf target, producing one `.trees` file per population.
   - For every population, chains three more gwf targets: `newick_simplify_vcf_remove` (tree/VCF extraction), `new_treva_smm` (likelihood inference), and `add_truth_values` (stamps the known true mutation rate and sample size onto the result).
   - Merges every population's result into a single `sim_likelihoods_smm.tsv`.

2. **`paralel.slim`**
   The SLiM simulation engine. If no input tree is given it founds a new ancestral population and seeds it with microsatellite mutations at random genomic positions; otherwise it loads the parent population's tree sequence, splits off a new descendant subpopulation of the requested size, and continues forward. A `modifyChild()` callback mutates each individual's microsatellites every generation, drawing the new repeat count from a transition matrix (`transition_matrix.csv`) built by `transitionMatrixMSprimeMethod.py` at the repo root. At the end of the run it saves the tree sequence (with STR metadata) so a descendant population can pick up from it.

3. **`modify_trees_get_vcf_infer_arg.py`**
   Post-processes one population's simulated tree sequence: recapitates and simplifies it down to the requested sample size, writes a per-locus Newick tree plus a (raw and allele-sorted) STR genotype table for every microsatellite, strips the STR mutations back out of the tree sequence, and writes the resulting clean tree sequence out as a VCF for downstream ARG inference.

4. **`new_treva_SMM.R`**
   The TREVA inference step. Reads the Newick trees and STR genotype table produced above and fits both an SMM (1 free parameter: mutation rate) and a TPM (3 free parameters: mutation rate, `p`, `m`) model via simulated annealing (`GenSA`), either per tree (`local`) or jointly across all trees (`global`). For each tree it keeps whichever model has the lower AIC and writes the result to `Overall_best_likelihood.tsv`.

## 2. Dependencies

- [gwf](https://gwf.app/) and [SLiM](https://messerlab.org/slim/) — see the repo-root README for installation and cluster backend setup.
- Python: `tskit`, `msprime`, `pyslim`, `pandas`, `numpy`
- R: `phytools`, `expm`, `ape`, `dplyr`, `stringr`, `GenSA`

## 3. Input

`workflow.py` expects a tab-separated metadata table (edge, parent, population size `N`, generations `gens`) describing the demographic tree — one row per population, with a `parent` column linking each population to the one it split from.

## 4. Running

```bash
gwf status
gwf run
```

This readme is the only AI generated information in this repository. All other files are written by the authors of this repository.. I'm sorry I was lazy. I'll try to updated.
