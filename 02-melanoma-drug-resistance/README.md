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
6. **DE pseudobulk (revisão estatística)** — o DE do passo 5a compara células,
   mas o tratamento é uma característica do tumor; as células não são réplicas
   independentes. Este passo soma as contagens brutas das células malignas por
   amostra (amostras com ≥ 10 células malignas: 12 treatment.naive, 8
   post.treatment) e roda DESeq2 com desenho `~ Cohort + treatment.group`
   (`06_pseudobulk_de.R`)

## 5. Results

- **Clustering:** 29 clusters não supervisionados (resolução 0.8) sobre
  7.186 células — granularidade bem mais fina que as ~7-8 categorias de
  `cell.types` publicadas (ver `cluster_vs_celltype_crosstab.csv`);
  resolução mais baixa (~0.4-0.5) pode aproximar melhor dos tipos
  celulares publicados, não testado nesta rodada.
- **DE em células malignas, post.treatment vs treatment.naive, nível de
  célula** (`de_malignant_post_vs_naive.csv`, script 05): **MIA** (avg log2FC
  = -3.65, p_adj = 4.7e-138) e **EMP1** (avg log2FC = -4.58, p_adj = 2.2e-135)
  caem pós-tratamento; **SERPINF1** (avg log2FC = +2.76) e **PFN1** (avg
  log2FC = +1.47) sobem. 5.384 genes passam p_adj < 0.05.
- **Mesmo contraste no nível de amostra (pseudobulk, script 06)**
  (`results/pseudobulk_de_malignant_post_vs_naive.csv`): apenas **117 genes**
  passam padj < 0.05, e só 61 dos 5.384 achados célula a célula (1,1%) se
  confirmam. Dos quatro genes destacados acima, **só EMP1 se mantém**
  (log2FC = -4.51, padj = 2.9e-06). MIA (padj = 0.11), SERPINF1 (padj = 0.32)
  e PFN1 (padj = 0.67) mantêm a direção mas não sobrevivem à correção com 8 vs
  12 tumores (`results/pseudobulk_vs_cell_highlight_genes.csv`,
  `figures/pseudobulk_highlight_genes.png`).
- **Score de resistência x fração de células T por amostra**
  (`resistance_score_vs_tcell_fraction.csv`): correlação de Spearman
  rho = -0.147, p = 0.50 (não significativa) — direção consistente com a
  hipótese (score mais alto → menos infiltração de T), mas sem poder
  estatístico suficiente com n = 31 amostras.

## 6. Interpretation

Os p-valores da análise célula a célula (script 05) estão inflados por
pseudo-replicação: tratam milhares de células de poucos tumores como
observações independentes. A análise pseudobulk (script 06) é a que deve ser
usada para inferência. Nela, a queda de **EMP1** pós-tratamento é robusta; a
queda de **MIA** e o aumento de **SERPINF1** aparecem como tendência na mesma
direção, mas sem significância após correção. O padrão é compatível com
desdiferenciação (programa MITF-low/AXL-high descrito na literatura de
resistência em melanoma), mas com estes dados a evidência sustenta EMP1, não o
programa completo. Isso está alinhado com a hipótese central
do paper original (Jerby-Arnon et al.) de um programa maligno intrínseco
ligado a evasão imune, embora aqui a análise seja simplificada (DE por
timing de biópsia, não o score de programa resistente completo do paper).

A correlação não significativa entre o score de resistência e a fração de
células T por amostra não permite concluir associação com exclusão imune
neste recorte dos dados — mais provável limitação de poder estatístico
(n=31) do que ausência real do efeito, já que a direção do efeito bate
com a hipótese.

## 7. Decision

**EMP1** é o candidato mais sólido (sobrevive à análise no nível do tumor);
MIA e SERPINF1 ficam como hipóteses secundárias, e PFN1 sai da lista. Os
candidatos são hipóteses razoáveis para priorizar em validação funcional (ex.: silenciar
MIA/EMP1 em linhagem de melanoma e avaliar reconhecimento por células T
in vitro) antes de qualquer decisão de desenvolvimento. Não há evidência
aqui forte o bastante para decisão de investimento em alvo — o achado
principal desta rodada é a viabilidade do pipeline e um sinal
qualitativamente coerente com a literatura, não uma confirmação
estatística robusta.

## 8. Limitations

- Não há classificação de resposta clínica (R/NR) neste dataset — a
  análise usa proxies (timing de tratamento, composição celular por
  amostra), não outcome clínico direto
- Mistura duas coortes (`Tirosh`, `New`) com possíveis diferenças
  técnicas/batch entre elas — não corrigido por batch integration nesta
  versão inicial
- N pequeno de tumores (31); correlação amostra-nível entre score de
  malignant program e fração de células T tem poder estatístico limitado
- O DE célula a célula (script 05) sofre de pseudo-replicação; usar o
  pseudobulk (script 06) para inferência. Com 8 vs 12 tumores o poder é
  limitado, e genes com log2FC muito alto no pseudobulk (ex.: MAGEA4) podem
  refletir poucos tumores com expressão extrema; inspecionar por amostra antes
  de interpretar
- O score de resistência do passo 5b foi construído com genes do DE célula a
  célula; refazer com os genes do pseudobulk antes de interpretar a correlação
  com a fração de células T
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
  `05_resistance_program_analysis.R` → `06_pseudobulk_de.R`
- Dados brutos/processados não versionados no git; regenerados a partir
  do script 01
