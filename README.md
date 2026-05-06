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

To isolate effects mediated by the **5-HT2A / IP3R2 signaling pathway**, results focus on features consistently changed in **both AV and IP relative to WT**. AV blocks the 5-HT2A receptor; IP lacks IP3R2-mediated calcium release. Effects shared between both groups implicate specific pathway. All values are Δ FC vs WT (*** p<0.001, ** p<0.01, * p<0.05, ns = not significant, BH FDR).

### Drug condition

| Feature | AV Q50 | AV Q75 | AV Q90 | IP Q50 | IP Q75 | IP Q90 |
|---------|--------|--------|--------|--------|--------|--------|
| Max simultaneous events | −0.91 *** | −2.13 *** | −6.41 *** | −0.95 *** | −2.23 *** | −6.70 *** |
| Area | −0.24 *** | −0.70 *** | −1.02 ** | ns | −1.16 *** | −2.32 *** |
| Perimeter | −0.16 ** | −0.41 *** | −0.55 *** | ns | −0.61 ** | −1.04 *** |
| Co-location (same location) | −0.33 *** | −0.33 ** | −0.33 ** | −0.33 * | ns | −0.67 ** |
| Co-location (similar size) | ns | −0.50 *** | −2.00 *** | ns | ns | −2.00 ** |

**Max simultaneous events** is the strongest finding. The suppression grows dramatically from Q50 (~−0.93 in both groups) to Q90 (~−6.5), meaning psilocybin generates extreme high-synchrony bursting in WT that is almost entirely absent in both AV and IP. This indicates complete dependence of synchronized network bursting on 5-HT2A/IP3R2 signaling.

**Area and perimeter** suppression is concentrated at Q75–Q90. Importantly, the IP effect at Q50 is not significant — smaller events persist without IP3R2 signaling, but large events are selectively absent. This indicates IP3R2 specifically drives high-amplitude calcium events rather than uniformly shifting event size.

**Co-localized events** show a tail-specific pattern: baseline-level spatial clustering (Q50) is intact in both groups, but high co-location counts at Q75–Q90 are strongly suppressed. Only the extreme clustering events depend on 5-HT2A/IP3R2 signaling.

**AV-only effects (not shared with IP):** Rise time, decay time, and duration 50–50% were significantly prolonged in AV at Q50 (all p<0.01) but showed no significant IP effects at any quantile. This may reflect additional 5-HT2A-mediated kinetic effects beyond IP3R2, or limited power in the IP group (~15 events/slice).

**No significant effects in either group:** Max ΔF, Max ΔF/F, all three AUC variants (raw/ΔF/ΔF/F), and duration 10–10% showed no significant effects after FDR correction. Psilocybin via 5-HT2A/IP3R2 selectively targets event size, spatial clustering, and network synchrony rather than uniformly altering all calcium event properties.

### Washout condition

| Feature | AV Q50 | AV Q75 | AV Q90 | IP Q50 | IP Q75 | IP Q90 |
|---------|--------|--------|--------|--------|--------|--------|
| Max simultaneous events | −1.13 *** | −3.54 *** | −7.26 *** | −1.46 *** | −4.00 *** | −7.80 *** |
| Area | −0.16 * | −0.45 ** | −0.80 ** | ns | ns | ns |
| Perimeter | ns | −0.32 *** | −0.39 ** | ns | ns | ns |
| Co-location (same location) | −0.25 *** | ns | −0.50 *** | ns | ns | ns |
| Co-location (similar size) | ns | −0.50 *** | −2.00 *** | ns | −1.00 *** | −2.00 ** |

Because washout values share the same baseline denominator as drug values, absolute Q50 values directly show recovery status (y = 1.0 = fully recovered):

| Feature | WT washout | AV drug | AV washout | IP drug | IP washout |
|---------|------------|---------|------------|---------|------------|
| Max simultaneous events | 2.13 | 2.13 | 1.00 | 1.18 | 0.67 |
| Area | 1.05 | 0.93 | 0.89 | 0.90 | 0.85 |
| Co-location (similar size) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

**None of the shared effects reverse at washout.** WT itself remains elevated above baseline for synchrony (Q50 = 2.13), showing psilocybin's effect persists in the intact pathway. AV recovers to baseline at Q50 for synchrony (washout Q50 = 1.00) but maintains strong suppression at Q75–Q90. Area remains partially suppressed in AV (washout Q50 = 0.89, below baseline) but IP loses significance at washout, likely due to reduced power (n = 58 washout events).

**IP shows progressive suppression rather than recovery.** Max simultaneous events Q50 drops from 1.18 during drug to 0.67 at washout — below baseline — while AV moves toward 1.0. IP3R2 knockout astrocytes do not reverse the network synchrony effect of psilocybin during the washout window, suggesting sustained downstream changes beyond acute receptor occupancy.

**Co-location (similar size)** maintains its tail-specific pattern at washout unchanged: Q50 = 1.00 in both groups (at baseline), Q75–Q90 strongly suppressed in both. The upper tail of spatial clustering does not recover.

**New washout-specific changes:** IP showed a rebound in Max ΔF/F (elevated above WT at all quantiles, p<0.001) not present during drug, suggesting compensatory upregulation of calcium amplitude. AV showed elevated Max ΔF above WT at Q50–Q90 (p≤0.023), also absent during drug — a rebound overshoot after receptor blockade is removed.

---

## Exported Files

All exports are saved under `Output__/` after running `analysis_revised.ipynb`.

**Per-group feature CSVs** — saved under `Output__/{group_path}/` for each group:

| File | Description |
|------|-------------|
| `{group}_baseline_normalized.csv` | Fold-change normalized, baseline condition |
| `{group}_drug_normalized.csv` | Fold-change normalized, drug condition |
| `{group}_washout_normalized.csv` | Fold-change normalized, washout condition |
| `{group}_baseline_unnormalized.csv` | Raw values, baseline |
| `{group}_drug_unnormalized.csv` | Raw values, drug |
| `{group}_washout_unnormalized.csv` | Raw values, washout |

**Per-group binned CSVs** — same location, one file per condition × bin size (9 files per group):

| File | Description |
|------|-------------|
| `{group}_{condition}_binned5.csv` | Mean per 5-frame bin |
| `{group}_{condition}_binned10.csv` | Mean per 10-frame bin |
| `{group}_{condition}_binned20.csv` | Mean per 20-frame bin |

**Event count exports** — saved directly under `Output__/`:

| File | Description |
|------|-------------|
| `all_groups_event_count_summary.csv` | One row per slice: raw counts + FC for all groups and conditions |
| `{group}_event_counts_raw.csv` | Frame-wise event counts per slice |
| `{group}_event_counts_raw_binned5.csv` | Frame counts summed per 5-frame bin |
| `{group}_event_counts_raw_binned10.csv` | Frame counts summed per 10-frame bin |
| `{group}_event_counts_raw_binned20.csv` | Frame counts summed per 20-frame bin |

Total: 12 feature CSVs + 9 binned CSVs per group, plus 17 event count files across all groups.
