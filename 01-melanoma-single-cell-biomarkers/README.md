# Projeto 01 — Single-Cell Biomarker Discovery in Melanoma

## 1. Problem

Nem todos os pacientes com melanoma metastático respondem a checkpoint
inhibitors (anti-PD-1/anti-CTLA-4). Identificar, a partir de dados de
scRNA-seq do microambiente tumoral, populações celulares e assinaturas
transcricionais associadas a resposta (responder, R) vs. não-resposta
(non-responder, NR) à imunoterapia.

## 2. Why it matters

Biomarcadores preditivos de resposta permitem estratificar pacientes antes
do tratamento, evitando exposição a toxicidade de imunoterapia sem benefício
e orientando combinações terapêuticas racionais (ex.: alvos em células
exaustas/Treg enriquecidas em NR). Isso é o núcleo de decisões de
desenvolvimento clínico e de medicina de precisão em oncologia.

## 3. Dataset

- **Accession:** [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978)
- **Referência:** Sade-Feldman et al., *Cell* 2018 — "Defining T Cell States
  Associated with Response to Checkpoint Immunotherapy in Melanoma"
- **Plataforma:** Smart-seq2 (full-length scRNA-seq, não é 10x/UMI-based)
- **Tamanho:** ~16.291 células de 48 biópsias tumorais (32 pacientes,
  amostras pré- e pós-tratamento)
- **Arquivos-chave na GEO:**
  - `GSE115978_counts.csv.gz` — matriz de contagens brutas
  - `GSE115978_tpm.csv.gz` — matriz normalizada (TPM)
  - `GSE115978_cell.annotations.csv.gz` — anotação de tipo celular e
    resposta clínica (R/NR) por célula, fornecida pelos autores

## 4. Methods

1. **Download** dos arquivos supplementary da GEO (`01_download_gse115978.R`)
2. **QC** — filtragem por número de genes detectados/célula, contagem total,
   % mitocondrial; remoção de outliers (`02_qc.R`)
3. **Pré-processamento Seurat** — criação do objeto, normalização
   (log-transform, já que os dados de origem são TPM/full-length e não
   UMI counts), seleção de genes variáveis, scaling, PCA
   (`03_seurat_preprocessing.R`)
4. **Clustering e anotação** — clustering (Louvain/Leiden via Seurat),
   UMAP, validação cruzada com as anotações de tipo celular publicadas
   (`04_clustering_annotation.R`)
5. **Análise de biomarcadores** — expressão diferencial R vs. NR por
   população celular (Wilcoxon/MAST), volcano plots, candidatos a
   assinatura preditiva (`05_biomarker_analysis.R`)

## 5. Results

*(a preencher após execução da análise)*

## 6. Interpretation

*(a preencher — interpretação biológica dos achados: quais populações e
genes se destacam, plausibilidade mecanística)*

## 7. Decision

*(a preencher — que decisão de desenvolvimento/pesquisa esses resultados
sustentariam, ex.: priorizar alvo X para validação funcional)*

## 8. Limitations

- Estudo observacional, n pequeno de pacientes (32); confundido por
  variáveis clínicas não totalmente controladas (linha de tratamento,
  local da biópsia)
- Smart-seq2 captura menos células por amostra que plataformas 10x,
  limitando poder estatístico por subtipo raro
- Anotações de resposta clínica (RECIST) podem não refletir perfeitamente
  a biologia tumoral no momento da biópsia
- Achados são hipótese-geradores; não validados experimentalmente aqui

## 9. Next steps

- Validar biomarcadores candidatos em coorte independente (ex.: dados de
  outro estudo de scRNA-seq em melanoma sob imunoterapia)
- Explorar deconvolução em dados bulk RNA-seq de coortes maiores usando a
  assinatura derivada
- Ligação com o Projeto 02 (drug resistance) para investigar sobreposição
  entre programas de resistência e não-resposta a imunoterapia

## 10. Reproducibility

- Ambiente: R (RStudio local) — ver [requirements.txt](requirements.txt)
  para versões de pacotes
- Seed fixa (`set.seed(42)`) em todos os scripts com componente estocástico
- Ordem de execução: `analysis/01_download_gse115978.R` →
  `02_qc.R` → `03_seurat_preprocessing.R` → `04_clustering_annotation.R` →
  `05_biomarker_analysis.R`
- Dados brutos/processados não versionados no git (ver `.gitignore`);
  regenerados a partir do script 01
