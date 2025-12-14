Protein–ligand docking을 위한 구조 데이터 전처리 
(RCSB PDB API 기반 구조 선별 및 PubChem API 기반 ligand curation 자동화)

## Overview
This project implements a preprocessing pipeline designed to curate
**docking-ready protein–ligand structures** prior to large-scale docking experiments.
Rather than focusing on docking itself, the pipeline addresses **input data quality**,
which is one of the most common causes of docking failure and poor reproducibility.

---

## Motivation
In large-scale docking workflows, failures frequently arise from:
- Low-resolution or inappropriate PDB structures
- Inconsistent ligand annotations across PDB entries
- Ambiguous protein–ligand relationships

To mitigate these issues, this pipeline was built to **systematically filter,
standardize, and organize protein–ligand structures** before docking is performed.

---

## What This Pipeline Does
- Automated retrieval of protein structures using the **RCSB PDB API**
- Selection of **X-ray crystal structures** with resolution-based quality filtering
- Parsing of ligand reports to identify **dockable ligands**
- Ligand curation using **PubChem API** to ensure consistency
- Construction of standardized **protein–ligand pairs** suitable for batch docking

This preprocessing step produces a curated dataset that can be directly used
for downstream docking and MD simulations.
