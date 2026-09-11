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

- **Accession:** [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575)
- **Referência:** Sade-Feldman et al., *Cell* 2018 — "Defining T Cell States
  Associated with Response to Checkpoint Immunotherapy in Melanoma"
  (PMID 30388456, DOI [10.1016/j.cell.2018.10.038](https://doi.org/10.1016/j.cell.2018.10.038))
- **Plataforma:** Smart-seq2 (full-length scRNA-seq), células CD45+ triadas
  por FACS
- **Tamanho:** 16.291 células de 48 biópsias tumorais (32 pacientes,
  amostras pré- e pós-tratamento com anti-CTLA4/anti-PD1/combinação)
- **Arquivos-chave na GEO (supplementary da série):**
  - `GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz` — matriz
    de expressão (genes x células), valores em log2[(TPM/10)+1]; **as 2
    primeiras linhas são cabeçalho** (linha 1 = ID da célula, linha 2 =
    ID do paciente/timepoint, ex. `Pre_P1`) — não são genes
  - `GSE120575_patient_ID_single_cells.txt.gz` — arquivo no formato de
    template de submissão GEO (não é um CSV simples); contém, a partir da
    linha 20, uma linha por célula com colunas: `Sample name`, `title`
    (= ID da célula, bate com o header da matriz TPM), `source name`,
    `organism`, `characteristics: patient ID (Pre=baseline; Post=on
    treatment)`, `characteristics: response` (**Responder**/**Non-responder**),
    `characteristics: therapy` (anti-CTLA4/anti-PD1/anti-CTLA4+PD1)
- **Nota de reprodutibilidade:** o rótulo de resposta clínica (R/NR) está
  neste arquivo por célula — **não** confundir com o GSE115978 usado no
  Projeto 02, que não tem esse rótulo (ver README do Projeto 02 para o
  histórico dessa confusão de accessions).

## 4. Methods

1. **Download** dos arquivos supplementary da GEO (`01_download_gse120575.R`)
2. **Parsing de metadados** — extração do bloco por célula do arquivo de
   template GEO (`patient_ID_single_cells.txt`), obtendo `cell_id`,
   `patient_id`, `response`, `therapy` (`02_parse_metadata.R`)
3. **QC** — filtragem por número de genes detectados/célula (TPM > 0) e
   % mitocondrial; a matriz já vem log-transformada, então os thresholds
   são adaptados a essa escala (`03_qc.R`)
4. **Pré-processamento Seurat** — criação do objeto a partir da matriz já
   log2-normalizada (sem reaplicar `NormalizeData` log), seleção de genes
   variáveis, scaling, PCA (`04_seurat_preprocessing.R`)
5. **Clustering e biomarcadores** — clustering (Louvain via Seurat), UMAP,
   e expressão diferencial Responder vs. Non-responder global e por
   cluster (`05_clustering_and_biomarkers.R`)

## 5. Results

*(a preencher após execução da análise)*

## 6. Interpretation

*(a preencher)*

## 7. Decision

*(a preencher)*

## 8. Limitations

- Estudo observacional, n pequeno de pacientes (32 pacientes / 48
  biópsias); confundido por variáveis clínicas (linha de terapia:
  anti-CTLA4 vs anti-PD1 vs combinação — não estratificado na análise
  inicial)
- Smart-seq2 captura menos células por amostra que plataformas 10x,
  limitando poder estatístico por subtipo raro
- Células são CD45+ apenas (compartimento imune) — não inclui células
  malignas/estromais (diferente do GSE115978 usado no Projeto 02)
- Matriz de expressão já vem log-transformada e sem contagens brutas
  disponíveis publicamente — métricas de QC (nº de genes, sem nº de
  reads) são mais limitadas que com dados de contagem bruta
- Achados são hipótese-geradores; não validados experimentalmente aqui

## 9. Next steps

- Cruzar biomarcadores candidatos com o programa de resistência de
  células malignas identificado no Projeto 02 (overlap entre não-resposta
  e exclusão de células T)
- Validar candidatos usando o marcador TCF7 destacado no paper original
  como preditor em coorte independente
- Estratificar a análise por linha de terapia (anti-CTLA4 vs anti-PD1 vs
  combinação), já que a resposta a cada uma pode ter biologia distinta

## 10. Reproducibility

- Ambiente: R (RStudio local) — ver [requirements.txt](requirements.txt)
  para versões de pacotes
- Seed fixa (`set.seed(42)`) em todos os scripts com componente estocástico
- Ordem de execução: `analysis/01_download_gse120575.R` →
  `02_parse_metadata.R` → `03_qc.R` → `04_seurat_preprocessing.R` →
  `05_clustering_and_biomarkers.R`
- Dados brutos/processados não versionados no git (ver `.gitignore`);
  regenerados a partir do script 01
- **Aviso de tamanho:** a matriz TPM descompactada tem ~4,5 GB — garanta
  espaço em disco livre antes de rodar o script 01
