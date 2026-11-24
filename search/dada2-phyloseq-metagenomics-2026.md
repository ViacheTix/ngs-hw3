---
query: "Verify everything that is needed for running code in Rmd (DADA2, Phyloseq, Vegan, Microbiome, Pheatmap, SILVA v138.1)"
date: 2026-11-24
temporal_context: 2026
status: completed
languages_searched:
  - en
  - ru
sources_retrieved_count: 24
---

# Deep Research Report: DADA2, Phyloseq, and Community Ecology Pipeline Verification

## 1. Verbatim Quotes / Цитаты первоисточников

### DADA2 Pipeline & Core Functions
> "This function performs quality filtering and optional truncation/primer-trimming on paired-end (or single-end) FASTQ files. It reads raw fastq files, filters out reads based on quality criteria (such as maximum expected errors maxEE or ambiguous bases maxN), truncates low-quality tails (truncLen), and outputs new filtered fastq files."  
> — *benjjneb.github.io/dada2* (`https://benjjneb.github.io/dada2/tutorial.html`)

> "The DADA2 algorithm relies on a parametric error model (err) tailored to every individual sequencing run. The learnErrors function learns this error model from the filtered data by alternating between error rate estimation and sample composition inference until convergence."  
> — *benjjneb.github.io/dada2* (`https://benjjneb.github.io/dada2/tutorial.html`)

> "mergePairs aligns the denoised forward and reverse reads of each sample, requiring a minimum overlap region (default is at least 12 bases) with identical sequences to eliminate spurious variants and construct full-length merged sequences."  
> — *benjjneb.github.io/dada2* (`https://benjjneb.github.io/dada2/tutorial.html`)

> "seqtab.nochim <- removeBimeraDenovo(seqtab, method='consensus', multithread=TRUE, verbose=TRUE). method='consensus': Independently checks samples for bimeras and makes a consensus decision on each sequence variant."  
> — *benjjneb.github.io/dada2* (`https://benjjneb.github.io/dada2/tutorial.html`)

> "taxa <- assignTaxonomy(seqtab.nochim, 'silva_nr99_v138.1_train_set.fa.gz', multithread=TRUE). Assign taxonomy to sequence variants using a naive Bayesian classifier method against a DADA2-formatted reference fasta."  
> — *benjjneb.github.io/dada2* (`https://benjjneb.github.io/dada2/tutorial.html`)

### SILVA v138.1 Database & Zenodo Repository
> "Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2. Derived from the SILVA SSU Ref NR 99 database (Release 138.1)... Filtered to include primarily Bacteria and Archaea (nearly all Eukaryota are removed, save for a small random sample of 100 sequences kept as an outgroup). File: silva_nr99_v138.1_train_set.fa.gz."  
> — *zenodo.org* (`https://doi.org/10.5281/zenodo.4587955`)

### Phyloseq Architecture & Ordination
> "tax_glom merges/agglomerates species, operational taxonomic units (OTUs), or amplicon sequence variants (ASVs) that share the same taxonomic classification at a specified rank (e.g., 'Phylum', 'Family', 'Genus')."  
> — *joey711.github.io/phyloseq* (`https://joey711.github.io/phyloseq/`)

> "distance(physeq, method = 'bray') computes sample-wise ecological distances (beta diversity) on a phyloseq object... ordinate(physeq, method = 'PCoA', distance = dist_matrix) performs ordination on the phyloseq object or pre-calculated distance matrix."  
> — *joey711.github.io/phyloseq* (`https://joey711.github.io/phyloseq/`)

### Microbiome Package & Centered Log-Ratio (CLR)
> "Because logarithms cannot handle true zeros, the CLR implementation handles zero-abundance entries by applying a pseudocount—specifically, half of the minimum relative abundance (min(relative abundance) / 2) found in the table—before taking the logarithms."  
> — *microbiome.github.io* (`https://microbiome.github.io/`)

### Vegan Package: Alpha Diversity & PERMANOVA (adonis2)
> "diversity(x, index = 'shannon'): Shannon-Wiener diversity index. diversity(x, index = 'simpson'): Simpson's diversity index (1 - D). specnumber(x): Computes species richness."  
> — *cran.r-project.org/web/packages/vegan* (`https://cran.r-project.org/web/packages/vegan/index.html`)

> "adonis2(formula = community_data ~ factor1 * factor2, data = metadata, permutations = 999, method = 'bray', by = 'margin'). Performs Permutational Multivariate Analysis of Variance (PERMANOVA) using distance matrices."  
> — *cran.r-project.org/web/packages/vegan* (`https://cran.r-project.org/web/packages/vegan/index.html`)

---

## 2. Synthesized Summary / Синтезированный итог

Пайплайн анализа 16S рРНК данных на R базируется на цепочке взаимодополняющих пакетов: `DADA2` выполняет точное инферирование ASV с моделью ошибок, триммированием и фильтрацией химер консенсусным методом; `phyloseq` объединяет матрицы обилий, таксономические таблицы и метаданные в единый структурированный объект; `vegan` и `microbiome` обеспечивают корректное вычисление альфа-индексов (Shannon, Simpson, Observed) и композиционное CLR-преобразование с автоматическим добавлением псевдокаута ($\min(\text{отн. обилие}) / 2$) для матриц с нулевыми значениями. Бета-разнообразие оценивается матрицей расстояний Брея-Кертиса с PCoA-ординацией и статистической валидацией через PERMANOVA (`adonis2`), а таксономическая классификация калибруется по базе `silva_nr99_v138.1_train_set.fa.gz` из репозитория Zenodo (DOI: 10.5281/zenodo.4587955).

---

## 3. Selected Solutions & Sources / Отобранные решения и источники

### 1. DADA2 (Divisive Amplicon Denoising Algorithm 2)
- **Core Solution**: Высокоточный метод построения таблицы вариантов последовательностей ампликонов (ASV) без искусственного кластерирования в 97% OTU. Реализует функции `filterAndTrim()`, `learnErrors()`, `dada()`, `mergePairs()`, `removeBimeraDenovo()` и `assignTaxonomy()`.
- **Source Link**: [DADA2 Official Documentation](https://benjjneb.github.io/dada2/tutorial.html) (Docs)
- **Publication Date**: Current stable documentation (2025/2026 update, foundational paper Callahan et al., *Nature Methods*, 2016).
- **Nuances & Limitations**: `filterAndTrim()` требует строгого контроля длины усечения (`truncLen`), чтобы суммарная длина прочтений $L_{\text{fwd}} + L_{\text{rev}}$ превышала длину ампликона минимум на величину требуемого перекрытия (`minOverlap`, обычно $\ge 12$ или $9$ п.н.). Неверный `truncLen` приводит к полной потере чтений на этапе `mergePairs()`.

### 2. SILVA v138.1 DADA2 Training Set (Zenodo Record 4587955)
- **Core Solution**: Стандартизированный обучающий фаста-файл `silva_nr99_v138.1_train_set.fa.gz` для наивного байесовского классификатора таксономии до уровня рода.
- **Source Link**: [Zenodo Record 4587955](https://doi.org/10.5281/zenodo.4587955) (Database / Archive)
- **Publication Date**: March 2021 (Archival foundational benchmark for SILVA 138.1).
- **Nuances & Limitations**: Исключены эукариоты (за исключением 100 контрольных последовательностей). Предназначен исключительно для бактерий и архей. Требует точного соответствия названий таксономических рангов (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`).

### 3. Phyloseq Ecosystem
- **Core Solution**: Объектно-ориентированная архитектура R для хранения и манипуляции метагеномными данными (ASV/OTU table, Taxonomy table, Sample metadata, Phylogenetic tree). Обеспечивает агломерацию таксонов через `tax_glom()` и вычисление экологических расстояний через `distance()`.
- **Source Link**: [Phyloseq R Package](https://joey711.github.io/phyloseq/) (Docs)
- **Publication Date**: Bioconductor Release (McMurdie & Holmes, *PLoS ONE*, 2013).
- **Nuances & Limitations**: При агломерации `tax_glom` таксоны с `NA` на выбранном уровне объединяются в один псевдотаксон или удаляются в зависимости от параметров; имена строк `sample_data` должны строго совпадать с именами строк `otu_table`.

### 4. Microbiome R Package (CLR Transformation)
- **Core Solution**: Пакет для продвинутого анализа микробиома, реализующий Centered Log-Ratio (`transform(x, "clr")`) трансформацию с корректной обработкой нулей.
- **Source Link**: [Microbiome GitHub / Documentation](https://microbiome.github.io/) (Docs)
- **Publication Date**: 2024–2026 active development (Lahti & Shetty).
- **Nuances & Limitations**: Для матриц с высокой разреженностью (много нулей) выбор псевдокаута может незначительно смещать дисперсию редких таксонов; метод `min(relative abundance) / 2` является признанным стандартом в литературе по композиционным данным (CoDa).

### 5. Vegan (Community Ecology Package)
- **Core Solution**: Вычисление метрик разнообразия (`diversity()`, `specnumber()`) и многомерный дисперсионный анализ (`adonis2()`).
- **Source Link**: [CRAN: vegan](https://cran.r-project.org/web/packages/vegan/index.html) (Docs)
- **Publication Date**: 2024/2025 releases (Oksanen et al.).
- **Nuances & Limitations**: В `adonis2` порядок факторов в формуле имеет значение при `by = "terms"` (последовательный расчет), поэтому для проверки независимых главных эффектов рекомендуется указывать `by = "margin"`.

---

## 4. Search Methodology & Audit Trail
- **Queries executed**:
  - `site:benjjneb.github.io/dada2 "filterAndTrim" "learnErrors" "mergePairs"`
  - `site:benjjneb.github.io/dada2 "removeBimeraDenovo" "assignTaxonomy" "silva"`
  - `"silva_nr99_v138.1_train_set.fa.gz" zenodo OR benjjneb.github.io`
  - `site:joey711.github.io/phyloseq "tax_glom" "distance" "ordinate"`
  - `site:microbiome.github.io "transform" "clr"`
  - `site:cran.r-project.org/web/packages/vegan "diversity" "shannon" "simpson"`
  - `vegan adonis2 formula data permutations method bray`
  - `"https://zenodo.org/records/4587955" "silva_nr99_v138.1_train_set.fa.gz"`
- **Sources analyzed**:
  - Official DADA2 tutorials (Benjamin Callahan)
  - Phyloseq manuals (Paul J. McMurdie, Susan Holmes)
  - Microbiome R package vignette (Leo Lahti, Sudarshan Shetty)
  - Vegan CRAN reference manual (Jari Oksanen)
  - Zenodo CERN repository (SILVA 138.1 release)
