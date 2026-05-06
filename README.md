# astrocyte-calcium-pipeline

Post-processing pipeline for astrocyte calcium imaging data exported from [AQuA2](https://github.com/yu-lab-vt/AQuA2). Handles per-slice normalization, quantile regression, and visualization across experimental groups.

---

## Notebooks

**`analysis.ipynb`** — run first
- Loads and transposes AQuA2 CSV exports
- Normalizes all features to baseline via per-slice fold-change
- Quantile regression at Q50 / Q75 / Q90 (WT as reference)
- Between-group dot plots
- Within-group timepoint plots (Baseline → Drug → Washout)
- Exports normalized CSVs and raw/binned event count CSVs

**`analysis_time.ipynb`** — run after main
- Loads normalized CSVs produced by the main notebook
- Bins events into 10-frame windows (frames 20–100)
- Events summarized within each slice per bin first, then aggregated across slices — giving each slice equal weight regardless of event count
- Time-series plots: group median ± MAD, all groups overlaid, Drug and Washout side by side
- Event count panel: FC-normalized per slice to its own baseline total

---

## Setup

**Dependencies:** `numpy pandas matplotlib scipy statsmodels`

Organize AQuA2 exports under `Output__/`:

```
Output__/
├── WT/
│   ├── data1/
│   │   ├── slice_baseline
│   │   ├── slice_psi
│   │   ├── slice_washout
│   │   └── ...
│   └── data2/
├── Antagonist- Volinanserin/
├── IP3R2 cKO/
└── CalEx/
```

Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_FRAME` | `20` | Start frame |
| `MAX_FRAME` | `100` | End frame |
| `CHANNEL_SUFFIX` | `'Ch1'` | Channel suffix in filenames; set `''` to omit |
| `CONTROL_GROUP` | `'WT'` | Reference group for quantile regression |
| `QUANTILES` | `[0.50, 0.75, 0.90]` | Quantiles tested |
| `GROUP_SIZE` | `10` | Frames per time bin (temporal notebook) |
| `MIN_SLICES_PER_BIN` | `3` | Bins below this are faded in temporal plots |

---

## Experimental Groups

| Group | Description |
|-------|-------------|
| WT | Wild Type + Psilocybin |
| AV | 5-HT2A antagonist (Volinanserin) + Psilocybin |
| IP | IP3R2 conditional knockout + Psilocybin |
| CE | CalEx + Psilocybin |

Three conditions per group: Baseline → Drug → Washout.
Slices: WT n=11, AV n=8, IP n=6. CE had no complete recording triplets and is excluded from statistics.

---

## Normalization

All features are normalized per slice using fold-change relative to that slice's own baseline median:

```
normalized_value = event_feature_value / slice_baseline_median
```

Normalization is computed per slice before pooling across slices, so events from preparations with different absolute activity levels are on a common scale. Both drug and washout values use the same baseline denominator — so **y = 1.0 is a shared reference throughout**: above 1.0 means elevated relative to baseline, below 1.0 means suppressed, and a washout value of 1.0 means full recovery.

---

## Statistical Methods

**Primary test: quantile regression at Q50 / Q75 / Q90**

For each feature, quantile regression is fitted across all events with group as a categorical predictor and WT as the reference:

```
feature_FC ~ C(Group)   [WT as reference]
```

The coefficient for each group at quantile q is the difference in the q-th percentile relative to WT. Quantile regression is used rather than a single summary test because key effects in calcium imaging data tend to be concentrated in the upper tail — large events, high-synchrony bursts — and would be missed by median-only comparisons.

BH FDR correction (α = 0.05) is applied within each quantile × comparison block.

**Note on statistical unit:** regression is applied at the event level. Events within the same slice are not fully independent; p-values reflect event-level sample size rather than slice-level replication (n = 6–11 slices per group). Coefficient magnitude is the primary effect-size estimate.

---

## Key Results

To isolate effects mediated by the **5-HT2A / IP3R2 signaling pathway**, results focus on features consistently changed in **both AV and IP relative to WT**. AV blocks the 5-HT2A receptor; IP lacks IP3R2-mediated calcium release. Effects shared between both groups implicate specific pathway. All values are Δ FC vs WT (*** p<0.001, ** p<0.01, * p<0.05, ns = not significant).

### Drug condition

| Feature | AV Q50 | AV Q75 | AV Q90 | IP Q50 | IP Q75 | IP Q90 |
|---------|--------|--------|--------|--------|--------|--------|
| Max simultaneous events | −0.91 *** | −2.13 *** | −6.41 *** | −0.95 *** | −2.23 *** | −6.70 *** |
| Area | −0.24 *** | −0.70 *** | −1.02 ** | ns | −1.16 *** | −2.32 *** |
| Perimeter | −0.16 ** | −0.41 *** | −0.55 *** | ns | −0.61 ** | −1.04 *** |
| Co-location (similar size) | ns | −0.50 *** | −2.00 *** | ns | ns | −2.00 ** |

**Max simultaneous events** is the strongest finding. The suppression grows dramatically from Q50 to Q90, meaning psilocybin generates extreme high-synchrony bursting in WT that is almost entirely absent in both AV and IP. This indicates complete dependence of synchronized network bursting on 5-HT2A/IP3R2 signaling.

**Area and perimeter** & **Co-localized events** show a tail-specific pattern, Q75–Q90 are strongly suppressed. 

**AV-only effects (not shared with IP):** Rise time, decay time, and duration 50–50% were significantly prolonged in AV at Q50 (all p<0.01) but showed no significant IP effects at any quantile. This may reflect additional 5-HT2A-mediated kinetic effects beyond IP3R2, or limited power in the IP group (~15 events/slice).

**No significant effects in either group:** Max ΔF, Max ΔF/F, all three AUC variants (raw/ΔF/ΔF/F), and duration 10–10% showed no significant effects after FDR correction. Psilocybin via 5-HT2A/IP3R2 selectively targets event size, spatial clustering, and network synchrony rather than uniformly altering all calcium event properties.

