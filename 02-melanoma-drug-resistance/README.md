# Projeto 02 — Malignant Cell Resistance Program & T Cell Exclusion in Melanoma

## 1. Problem

Identificar, a partir de scRNA-seq de 31 tumores de melanoma, um programa
transcricional intrínseco às células malignas associado à exclusão de
células T do microambiente tumoral — um mecanismo de resistência primária
à imunoterapia que independe de mutação/perda de antígeno.

## 2. Why it matters

Um programa de resistência detectável em células malignas *antes* do
início do tratamento permitiria identificar, prospectivamente, pacientes
com baixa probabilidade de resposta a checkpoint inhibitors — e aponta
alvos combináveis (ex.: inibição de CDK4/6, testada no paper original)
para reverter a exclusão imune e resgatar a resposta à imunoterapia.

## 3. Dataset

- **Accession:** [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978)
- **Referência:** Jerby-Arnon et al., *Cell* 2018 — "Single-cell RNA-seq of
  melanoma ecosystems reveals sources of T cell exclusion linked to
  immunotherapy clinical outcomes" (PMID 30388455)
- **Plataforma:** Smart-seq2 (full-length scRNA-seq)
- **Tamanho:** ~7.186 células (após QC) de 31 tumores de melanoma;
  combina 2 coortes — `Tirosh` (dados de Tirosh et al. 2016, *Science*,
  GSE72056) e `New` (dados novos deste paper)
- **Arquivos-chave na GEO:**
  - `GSE115978_counts.csv.gz` — matriz de contagens brutas
  - `GSE115978_tpm.csv.gz` — matriz normalizada (TPM)
  - `GSE115978_cell.annotations.csv.gz` — anotação por célula: `cells`
    (ID), `samples`, `cell.types`, `treatment.group`
    (`post.treatment`/`treatment.naive`), `Cohort`, `no.of.genes`,
    `no.of.reads`

**Nota importante:** este dataset **não** contém rótulo de resposta
clínica (Responder/Non-responder) por célula ou por amostra — apenas
timing da biópsia (`treatment.group`). A classificação R/NR usada no
Projeto 01 vem de um dataset diferente (GSE120575, Sade-Feldman et al.).
Os dois papers foram publicados lado a lado na mesma edição da *Cell*
pelo mesmo grupo (Broad/MGH), o que gerou a confusão inicial de
accessions neste projeto — documentando aqui para não se repetir.

## 4. Methods

1. **Download** dos arquivos supplementary da GEO (`01_download_gse115978.R`)
2. **QC** — filtragem por genes detectados/célula, contagem total, % mito
   (`02_qc.R`)
3. **Pré-processamento Seurat** — normalização, genes variáveis, PCA
   (`03_seurat_preprocessing.R`)
4. **Clustering e anotação** — clustering, UMAP, validação cruzada com
   `cell.types` publicado (`04_clustering_annotation.R`)
5. **Programa de resistência x exclusão de células T** — em células
   malignas (`cell.types == "Mal"`), (a) DE entre `post.treatment` vs
   `treatment.naive` como proxy de programa induzido/selecionado por
   tratamento; (b) score de assinatura por amostra correlacionado com a
   fração de células T na mesma amostra (proxy de exclusão imune)
   (`05_resistance_program_analysis.R`)

## 5. Results

*(a preencher após execução da análise)*

## 6. Interpretation

*(a preencher)*

## 7. Decision

*(a preencher)*

## 8. Limitations

- Não há classificação de resposta clínica (R/NR) neste dataset — a
  análise usa proxies (timing de tratamento, composição celular por
  amostra), não outcome clínico direto
- Mistura duas coortes (`Tirosh`, `New`) com possíveis diferenças
  técnicas/batch entre elas — não corrigido por batch integration nesta
  versão inicial
- N pequeno de tumores (31); correlação amostra-nível entre score de
  malignant program e fração de células T tem poder estatístico limitado
- Direção causal não estabelecida — score alto de "programa de
  resistência" correlacionado com baixa fração de células T é consistente
  com, mas não prova, exclusão causada pelo programa maligno

## 9. Next steps

- Cruzar genes do programa de resistência identificado aqui com os
  candidatos a biomarcador do Projeto 01 (overlap entre "não-resposta" e
  "exclusão de células T")
- Validar o score de assinatura na coorte independente de 112 pacientes
  mencionada no paper original (dados bulk RNA-seq)
- Testar enriquecimento do programa para alvos de CDK4/6 e outras vias
  descritas no paper como moduladoras da exclusão imune

## 10. Reproducibility

- Ambiente: R (RStudio local) — ver [requirements.txt](requirements.txt)
- Seed fixa (`set.seed(42)`) em todos os scripts com componente estocástico
- Ordem de execução: `analysis/01_download_gse115978.R` →
  `02_qc.R` → `03_seurat_preprocessing.R` → `04_clustering_annotation.R` →
  `05_resistance_program_analysis.R`
- Dados brutos/processados não versionados no git; regenerados a partir
  do script 01
