# HTr_case_control

Statistical hit calling for HT-recruit (HTr) case/control screens.

Compares transcription activation domain activity between two conditions using a multinomial bootstrap approach. For each domain, read counts in the ON and OFF fractions are resampled to build a distribution of plausible log2(ON/OFF) enrichment values (rho). The difference in rho distributions between case and control conditions yields a posterior probability of dominance for each domain, calibrated against a null population of non-activating domains. Hits are called at user-specified false positive rates relative to the null.

---

## Installation

```bash
pip install git+https://github.com/lillianpetersen/HTr_case_control.git
```

Or in editable mode for development:

```bash
git clone https://github.com/lillianpetersen/HTr_case_control.git
cd HTr_case_control
pip install -e ".[dev]"
```

---

## Input data

Each condition (case and control) is a CSV file with one row per domain and the following count columns:

| Column | Description |
|---|---|
| `countsOFF_R1` | Reads in the OFF fraction, replicate 1 |
| `countsON_R1` | Reads in the ON fraction, replicate 1 |
| `countsOFF_R2` | Reads in the OFF fraction, replicate 2 |
| `countsON_R2` | Reads in the ON fraction, replicate 2 |

The row index must be a semicolon-delimited domain label. Two formats are parsed automatically:

- `mBR;library_type;protein_id;gene_name;tile_number;status`
- `TFTiles;ENSG_ID;gene_name;tile_index;start_aa;end_aa`

Parsed domain names take the form `GeneName;TileNumber`.

---

## Running the analysis

```bash
python run_package_test.py \
    --case    path/to/Treatment_combo.csv \
    --control path/to/Control_combo.csv \
    --output  path/to/output_folder/
```

Condition names are inferred from filenames by stripping `_combo.csv` (or `.csv`).

### Required arguments

| Argument | Description |
|---|---|
| `--case` | Path to the case condition CSV |
| `--control` | Path to the control condition CSV |
| `--output` | Output folder for all results CSVs and figures |

### Optional arguments

#### Null population

| Argument | Default | Description |
|---|---|---|
| `--null-file` | *(auto)* | CSV specifying the null population (see formats below). If omitted, null domains are classified automatically from the control condition. |

**Automatic null classification** (when `--null-file` is not provided):
1. Domains whose label contains `random_control` (case-insensitive) are always null.
2. Rho is computed for all domains in the control condition. Domains with `rho_mean ≤ mean(random_controls) + 2·SD` are classified as non-activating (null).
3. Domains above that threshold are activating and eligible to be called hits.

**`--null-file` formats:**
- **CRTF annotation**: columns `HGNC Symbol`, `Tile Number`, `Hit minCMV` — domains with `Hit minCMV == "Hit"` are activating.
- **is_null**: index = domain labels, column `is_null` (bool).
- **activating**: index = domain labels, column `activating` (bool).

#### Replicates

| Argument | Default | Description |
|---|---|---|
| `--drop-rep` | *(none)* | Drop replicate `1` or `2` from both conditions. Applied before analysis; both replicates are still shown on replicate scatter plots. |

#### Count filtering

| Argument | Default | Description |
|---|---|---|
| `--min-per-bucket` | `250` | Minimum reads in any single ON or OFF fraction per replicate. |
| `--min-per-replicate-total` | `700` | Minimum total ON+OFF reads per replicate. |

#### Bootstrapping

| Argument | Default | Description |
|---|---|---|
| `--n-boot` | `1000` | Bootstrap iterations per replicate. |
| `--add-counts` | `100` | Pseudocounts added per bootstrap draw, stabilising log-ratio estimates for low-count domains. |

#### Hit calling thresholds

Thresholds are set from quantiles of the null posterior probability distribution.

| Argument | Default | Description |
|---|---|---|
| `--fpr-strong` | `0.005` | Tail quantile for strong hits (0.5% FPR). |
| `--fpr-weak` | `0.02` | Tail quantile for weak hits (2% FPR). |

---

## Output

All files are written to `--output`.

### CSVs

| File | Description |
|---|---|
| `{Case}_vs_{Control}_results.csv` | One row per domain passing count filters. Contains rho scores, posterior probability, and hit call columns. |
| `{Case}_vs_{Control}_hits.csv` | Called hits only, sorted by confidence. |

**Key columns:**

| Column | Description |
|---|---|
| `{Case} R1/R2/Avg/SD` | Rho (log2 ON/OFF) per replicate, mean, and SD for the case condition |
| `{Control} R1/R2/Avg/SD` | Same for the control condition |
| `post_prob` | Posterior probability of dominance (0–1); >0.5 = more active in case, <0.5 = less active |
| `logit_dom` | logit(post_prob); used as the significance axis |
| `sig_up` / `sig_down` | Strong hit calls (bool) |
| `sig_up_weak` / `sig_down_weak` | Weak hit calls (bool) |
| `is_null` | True if domain belongs to the null population |

### Figures

| File | Description |
|---|---|
| `replicates_{Case}.pdf` | R1 vs R2 rho scatter for the case condition |
| `replicates_{Control}.pdf` | R1 vs R2 rho scatter for the control condition |
| `scatter_{Case}_vs_{Control}.pdf` | Case vs control rho coloured by hit call |
| `delta_plot_{Case}_vs_{Control}.pdf` | Δrho (case − control) vs baseline rho (control) |
| `delta_plot_{Case}_vs_{Control}_labeled.pdf` | Same with hit labels |
| `volcano_{Case}_vs_{Control}.pdf` | Δrho vs \|logit posterior\|; hits in upper corners |
| `ranked_strip_{Case}_vs_{Control}.pdf` | Domains ranked by logit posterior with top hits labelled |
| `logit_posterior_{Case}_vs_{Control}.pdf` | Histogram of logit posterior; null on secondary y-axis |
| `bootstrap_examples_{Case}.pdf` | Bootstrap delta distributions for example domains in each hit category |

---

## Examples

```bash
# Basic run
python run_package_test.py \
    --case    results/ERKi_combo.csv \
    --control results/Veh_combo.csv \
    --output  output/ERKi_vs_Veh/

# Drop a technically failed replicate
python run_package_test.py \
    --case    results/583TF_IFN_combo.csv \
    --control results/583TF_Veh_combo.csv \
    --output  output/583TF_IFN_vs_Veh/ \
    --drop-rep 2

# Use a CRTF annotation file for null classification
python run_package_test.py \
    --case      results/ERKi_combo.csv \
    --control   results/Veh_combo.csv \
    --output    output/ERKi_vs_Veh/ \
    --null-file annotations/CRTF_hits.csv

# Stricter hit calling with more bootstrap iterations
python run_package_test.py \
    --case       results/ERKi_combo.csv \
    --control    results/Veh_combo.csv \
    --output     output/ERKi_vs_Veh_strict/ \
    --n-boot     2000 \
    --fpr-strong 0.001 \
    --fpr-weak   0.01
```
