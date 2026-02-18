# ECG QRS-Based Beat Selection Pipeline (Script 02/03)

This repository contains an R script for QRS-based preprocessing of ECG feature data following time alignment to experimental events (e.g. drug injection).

This script represents the **second step (02/03)** of a modular ECG analysis pipeline.

It takes time-windowed ECG recordings as input and performs objective QRS cutoff detection followed by standardized First-N beat averaging, producing cleaned summary tables ready for statistical analysis.

---

## Overview

- Automatic QRS cutoff detection using a 2-component Gaussian mixture model
- Multiple cutoff strategies (Comp1 percentile, equal-density, valley, posterior probability)
- Interactive GUI (no command-line arguments required)
- Standardized First-N beat extraction after QRS threshold crossing
- Per-recording numeric summaries (mean values)
- Clean, analysis-ready outputs (CSV + XLSX)
- Robust to Windows / RStudio Tk issues

---

## What the pipeline does

1. Starting from per-recording ECG CSV files (e.g. post-injection 2-minute windows), the script:
2. Loads multiple ECG CSV files (one file per mouse / subject)
3. Pools QRS interval values across all recordings
4. Fits a 2-component Gaussian mixture model to the QRS distribution
5. Computes multiple candidate QRS cutoffs:
       - Comp1 percentiles
       - Equal-density intersection
       - Mixture valley
       - Posterior probability–based cutoff
6. Allows interactive selection of the cutoff strategy via GUI
7. For each recording:
       - Identifies the first beat exceeding the QRS cutoff
       - Extracts the next N beats (user-defined via textbox)
       - Computes the mean of all numeric ECG features
8. Exports cleaned summary tables with standardized column names
9. Writes all outputs to a structured, analysis-ready folder hierarchy

All file paths, cutoff strategies, and parameters are selected interactively via GUI dialogs.

---

## Typical use cases

- ECG drug-response analysis after time alignment
- Objective exclusion of early physiological beats
- Standardized beat selection across animals
- High-throughput ECG phenotyping
- Preparation of clean summary tables for figures
- Reproducible preprocessing in collaborative projects

---

## Position in the ECG pipeline

This script is Script 02 of a 3-step ECG analysis workflow:

1. Script 01 – Injection-aligned preprocessing and window generation
2. Script 02 (this repository) – QRS-based beat selection and First-N averaging
3. Script 03 – Metadata integration, statistics, and visualization

**Metadata integration (grouping, treatment, genotype, etc.) is intentionally performed outside this script and used in downstream analysis (Script 03).**

---

## Methods Description

ECG recordings were processed using a custom R pipeline in which QRS interval distributions were pooled across recordings and modeled using a two-component Gaussian mixture. An objective QRS cutoff was selected interactively, and for each recording the first beat exceeding this threshold was identified. The subsequent N beats were extracted and averaged to generate per-recording ECG feature summaries for downstream statistical analysis.

---

## Authorship

This script was developed by **Michele Buono** and can be used freely for research purposes, provided appropriate citation of the author.

**Buono, M. F. (2026). ECG Analysis Pipeline (Script 02/03) - QRS cutoff & First-N Averaging (v.1.0.0). Zenodo. 
https://doi.org/10.5281/zenodo.18504302**

The overall workflow, structure, and clarity of the pipeline were iteratively refined with assistance from **ChatGPT 5.2**, which was used as a tool to improve code organization, documentation, and usability.


