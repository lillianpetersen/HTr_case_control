# htr_hit_caller

Statistical hit calling for HT-recruit (HTr) screens.

---

## What is this?

**HT-recruit** is a high-throughput assay that measures whether a protein domain can activate transcription. A large library of candidate domains is expressed in cells, and sequencing tells you how many times each domain ended up in the **ON** fraction (activating) versus the **OFF** fraction (not activating). The ratio of ON to OFF counts — called **rho** — measures how strongly each domain activates transcription.

This package compares rho values between two experimental conditions — for example, cells treated with a kinase inhibitor versus cells treated with vehicle (DMSO) — to find domains whose activity **changed** due to the treatment. These are the **hits**.

### Why bootstrapping?

Sequencing counts are noisy. If you just compared average rho values between conditions, you would not know whether a difference is real or just random noise. This package uses **bootstrapping**: it resamples the count data thousands of times to build a distribution of plausible rho values for each domain. By comparing these distributions between conditions, it estimates a **posterior probability** for each domain — essentially, how confident we are that its activity truly changed.

- Posterior probability near **1** → strong evidence the domain became MORE active in the treatment condition
- Posterior probability near **0.5** → no clear change
- Posterior probability near **0** → strong evidence the domain became LESS active

---

## Input data

Each condition (case and control) is a CSV file with one row per domain and the following columns:

| Column | Description |
|---|---|
| `countsOFF_R1` | Read counts in the OFF fraction, replicate 1 |
| `countsON_R1` | Read counts in the ON fraction, replicate 1 |
| `countsOFF_R2` | Read counts in the OFF fraction, replicate 2 |
| `countsON_R2` | Read counts in the ON fraction, replicate 2 |
| `R1`, `R2`, `Avg`, `Standard Error` | Pre-computed summary columns (ignored) |

The row index (first column) must be a semicolon-separated label identifying each domain. Two label formats are supported automatically:

- **mBR library** (kinase inhibitor screen): `mBR;library_type;protein_id;gene_name;tile_number;status`
- **TFTiles library** (TF screen): `TFTiles;ENSG_ID;gene_name;tile_index;start_aa;end_aa`

---

## What is a "null domain"?

To call hits, the analysis needs to know which domains should NOT be hits — the **null population**. These are domains we expect to show no change between conditions. The null population serves as a baseline to calibrate the thresholds: if a domain's posterior probability is well outside the range seen in null domains, it is called a hit.

The null population is typically:
1. **Random sequence controls** — synthetic sequences not expected to activate transcription
2. **Non-activating domains** — domains that show very low activity in the control condition

---

## Installation

```bash
cd htr_hit_caller
pip install -e ".[dev]"
```

---

## Running the analysis

The main analysis script is `run_package_test.py` in the `LP012_kinase_inh_HTr_analysis` directory. Run it from the command line:

```bash
python run_package_test.py \
    --case    path/to/ERKi_combo.csv \
    --control path/to/Veh_combo.csv \
    --output  path/to/output_folder/
```

The condition names (e.g. `ERKi`, `Veh`) are automatically inferred from the filenames by stripping `_combo.csv`.

### Required arguments

| Argument | Description |
|---|---|
| `--case` | Path to the CSV file for the **treatment** condition (e.g. kinase inhibitor) |
| `--control` | Path to the CSV file for the **control** condition (e.g. vehicle/DMSO) |
| `--output` | Folder where all results files and plots will be saved |

### Optional arguments

#### Null population

| Argument | Default | Description |
|---|---|---|
| `--null-file` | *(auto)* | CSV specifying which domains belong to the null population (see below). If not provided, null domains are classified automatically from the control condition data. |

**Automatic null classification** (when `--null-file` is not provided):
1. Any domain whose label contains `random_control` is always null.
2. For all other domains, rho is computed from the control condition. Domains with rho below the mean + 2 standard deviations of the random controls are classified as non-activating (null).
3. The remaining domains — those above that threshold — are considered activating and are eligible to be called hits.

**Providing a null file** — three formats are accepted:
- **CRTF annotation format**: a CSV with columns `HGNC Symbol`, `Tile Number`, and `Hit minCMV`. Domains marked `Hit` in `Hit minCMV` are treated as activating; everything else is null.
- **is_null format**: a CSV whose index matches domain labels and has a column `is_null` (True = null).
- **activating format**: a CSV whose index matches domain labels and has a column `activating` (True = activating, False = null).

#### Replicates

| Argument | Default | Description |
|---|---|---|
| `--drop-rep` | *(none)* | Drop replicate `1` or `2` from **both** conditions. Use this when one replicate has a technical problem (e.g. compressed sequencing). The kept replicate is used for the full analysis; both replicates are still shown on the replicate scatter plots. |

#### Count filtering

Domains with very few sequencing reads are unreliable. These thresholds remove them before analysis.

| Argument | Default | Description |
|---|---|---|
| `--min-per-bucket` | `250` | Minimum number of reads in any single ON or OFF fraction. Domains below this in any replicate are dropped. |
| `--min-per-replicate-total` | `700` | Minimum total reads (ON + OFF combined) per replicate. Domains below this are dropped. |

If your library has lower sequencing depth than usual, you may need to lower these values. If you have very deep sequencing, you could raise them for stricter quality control.

#### Bootstrapping

| Argument | Default | Description |
|---|---|---|
| `--n-boot` | `1000` | Number of times the data is resampled per replicate. Higher values give more stable results but take longer to run. 1000 is usually sufficient. |
| `--add-counts` | `100` | Pseudocounts added to each bootstrap resample. These stabilise the log-ratio calculation for domains with low counts by pulling estimates slightly toward zero. |

#### Hit calling thresholds

Hits are called at two stringency levels. The thresholds are set so that the expected **false positive rate** (fraction of null domains incorrectly called as hits) equals the values below.

| Argument | Default | Description |
|---|---|---|
| `--fpr-strong` | `0.005` | False positive rate for **strong** hits (0.5%). Use stricter values (e.g. `0.001`) if you want fewer, higher-confidence hits. |
| `--fpr-weak` | `0.02` | False positive rate for **weak** hits (2%). These are suggestive but less certain. |

---

## Output files

All files are written to the folder specified by `--output`.

### Results CSVs

| File | Description |
|---|---|
| `{Case}_vs_{Control}_results.csv` | One row per domain passing count filters. Contains rho scores for both conditions, posterior probability, and hit call columns. |
| `{Case}_vs_{Control}_hits.csv` | Subset of results containing only called hits, sorted by confidence. |

**Key columns in the results file:**

| Column | Description |
|---|---|
| `{Case} R1`, `{Case} R2` | Rho (log2 ON/OFF) for each replicate in the case condition |
| `{Case} Avg`, `{Case} SD` | Mean and standard deviation of rho across replicates |
| `{Control} R1/R2/Avg/SD` | Same for the control condition |
| `post_prob` | Posterior probability of activity change (0–1). Values near 1 = more active in case; near 0 = less active in case; near 0.5 = no change. |
| `sig_up` | `True` if the domain is a **strong up** hit (more active in treatment) |
| `sig_down` | `True` if the domain is a **strong down** hit (less active in treatment) |
| `sig_up_weak` | `True` if the domain is a **weak up** hit |
| `sig_down_weak` | `True` if the domain is a **weak down** hit |
| `is_null` | `True` if the domain belongs to the null population |

### Plots

| File | Description |
|---|---|
| `replicates_{Case}.pdf` | Scatter plot of R1 vs R2 rho for the case condition. Points should fall near the diagonal if replicates agree well. |
| `replicates_{Control}.pdf` | Same for the control condition. |
| `scatter_{Case}_vs_{Control}.pdf` | Case vs control rho for every domain, coloured by hit call. Points above the diagonal became more active; below became less active. |
| `delta_plot_{Case}_vs_{Control}.pdf` | Change in activity (case − control) vs baseline activity (control rho). |
| `delta_plot_{Case}_vs_{Control}_labeled.pdf` | Same with hit labels. |
| `volcano_{Case}_vs_{Control}.pdf` | Volcano plot: x = change in activity, y = statistical confidence. Hits appear in the upper corners. |
| `ranked_strip_{Case}_vs_{Control}.pdf` | All domains ranked by confidence, with the top hits labelled. |
| `logit_posterior_{Case}_vs_{Control}.pdf` | Histogram of posterior probability scores. Null domains (right axis) should be spread evenly; hits cluster at the edges. |
| `bootstrap_examples_{Case}.pdf` | Bootstrapping distributions for a few example domains in each hit category. Illustrates how the analysis works. |

---

## Examples

**Basic run:**
```bash
python run_package_test.py \
    --case    results/ERKi_combo.csv \
    --control results/Veh_combo.csv \
    --output  output/ERKi_vs_Veh/
```

**Drop a bad replicate:**
```bash
python run_package_test.py \
    --case    results/583TF_IFN_combo.csv \
    --control results/583TF_Veh_combo.csv \
    --output  output/583TF_IFN_vs_Veh/ \
    --drop-rep 2
```

**Use a CRTF annotation file for null classification:**
```bash
python run_package_test.py \
    --case      results/ERKi_combo.csv \
    --control   results/Veh_combo.csv \
    --output    output/ERKi_vs_Veh/ \
    --null-file annotations/CRTF_hits.csv
```

**Stricter hit calling with more bootstrap iterations:**
```bash
python run_package_test.py \
    --case       results/ERKi_combo.csv \
    --control    results/Veh_combo.csv \
    --output     output/ERKi_vs_Veh_strict/ \
    --n-boot     2000 \
    --fpr-strong 0.001 \
    --fpr-weak   0.01
```

**Lower count filters for a lower-depth screen:**
```bash
python run_package_test.py \
    --case                   results/ERKi_combo.csv \
    --control                results/Veh_combo.csv \
    --output                 output/ERKi_vs_Veh/ \
    --min-per-bucket         100 \
    --min-per-replicate-total 400
```
