# Computational pharmacology portfolio

Analyses in translational oncology and computational pharmacology, built on public data.
Each folder is a self-contained project with its own README (question, data, methods,
results, limitations), scripts, figures and pinned package versions.

| Project | Question | Methods | Language |
|---|---|---|---|
| [melanoma-scrnaseq-resistance](melanoma-scrnaseq-resistance/) | Which malignant-cell genes change after immunotherapy in melanoma, once tumours rather than cells are treated as replicates? | Seurat, pseudobulk DESeq2 on 20 tumours (GSE115978) | R |
| [solubility-scaffold-shift](solubility-scaffold-shift/) | How much accuracy does a solubility model lose on chemical series it has never seen? | RDKit, scaffold split, random forest and XGBoost, bootstrap CIs (AqSolDB) | Python |
| [oncology-trials-brazil](oncology-trials-brazil/) | Which oncology trials run in Brazil, who sponsors them and where are the sites? | ClinicalTrials.gov API v2, pandas, state map with IBGE boundaries | Python |

Raw and processed data are not versioned; each project downloads them with its first
script. Results tables and figures are versioned so they can be read without rerunning
anything.

Related repositories: [scanpy-lowmem](https://github.com/thaisdolzany/scanpy-lowmem) (single-cell
analysis on a 16 GB laptop), [omics-fetch](https://github.com/thaisdolzany/omics-fetch)
(GEO, cBioPortal and TCGA downloads) and [ngs-germline-panel-nf](https://github.com/thaisdolzany/ngs-germline-panel-nf)
(Nextflow pipeline for germline variant calling in gene panels, benchmarked against a
truth set).

Author: Thaís Dolzany, PhD · [thaisdolzany.github.io](https://thaisdolzany.github.io) ·
[ORCID 0000-0002-5503-4051](https://orcid.org/0000-0002-5503-4051)

License: MIT
