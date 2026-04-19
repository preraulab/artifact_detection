# artifact_detection

Iterative z-score artifact detection for continuous time-series data (primarily sleep EEG). Produces a logical mask identifying samples that should be excluded from downstream analysis.

## Entry point

```matlab
artifacts = detect_artifacts(data, Fs, 'Name', Value, ...)
```

- `data`  — `Nx1 double`, single-channel time series (required)
- `Fs`    — scalar double, sampling frequency in Hz (required)
- Returns `1xN logical` — `true` at samples flagged as artifact

See the docstring in `detect_artifacts.m` for every name-value pair. The defaults are tuned for lab-standard sleep EEG at 100–500 Hz.

## Algorithm overview

The pipeline runs in four stages. The input is a 1-D time series; the output is a logical mask the same length as the input.

```
                    +---------------------------+
   raw data ------> |  Stage 1: Pre-filter      |
                    |   - NaN / Inf samples     |
                    |   - Flat stretches        |
                    |   - ±10σ outliers         |
                    |   - Spectral slope test   |
                    |   - Interpolate bad spots |
                    +------------+--------------+
                                 |
                     data_fixed (no NaN/gaps)
                                 |
            +--------------------+--------------------+
            |                                         |
            v                                         v
  +---------+---------+                    +----------+---------+
  | Stage 2: HF band  |                    | Stage 3: BB band   |
  | high-pass 35 Hz   |                    | high-pass 0.1 Hz   |
  | envelope -> log   |                    | envelope -> log    |
  | iterative z-score |                    | iterative z-score  |
  +---------+---------+                    +----------+---------+
            |                                         |
            +--------------------+--------------------+
                                 v
                    +---------------------------+
                    |  Stage 4: Union + buffer  |
                    |  hf | bb | pre-flagged    |
                    |  optional ±buffer secs    |
                    +------------+--------------+
                                 |
                             artifacts
```

---

### Stage 1 — Preliminary bad-index detection

Before the main detection loop runs, the input is screened for gross problems. Flagged samples become "seeds" that are (a) excluded from z-score statistics and (b) always marked artifact in the final output. The data is then gap-filled so the downstream filters see a continuous signal.

**1a. NaN / Inf samples** — flagged directly.

**1b. Flat regions** — `get_chunks(data, Fs)` returns, among other things, a logical mask `is_flat` marking stretches of identical consecutive values (constant signal, typically disconnected electrode or saturation). Exact threshold for "flat" is defined inside `get_chunks`.

**1c. Outlier samples** — after the above are excluded, compute `μ` and `σ` of the remaining good data. Flag any sample outside `μ ± 10σ`. Hard-coded multiplier (`outlier_scalar = 10`) handles gross amplitude glitches the iterative z-score might dilute.

**1d. Spectral slope test** *(enabled by default, `slope_test=true`)*

Healthy EEG has a negative log–log power slope (roughly 1/f). Broken channels, amplifier saturation, or electrical contamination flatten or reverse this slope.

- Compute a multitaper spectrogram of the whole recording:
  - Frequency range: `[1, min(55, Fs/2)]` Hz
  - Time–bandwidth product: 10, 19 tapers (maximum recommended for TW=10)
  - Window: 10 s, step 5 s
  - Uses `multitaper_spectrogram_mex` (compiled MEX for speed)
- For each time window, fit `log(power) = a + b·log(freq)` by `polyfit(..., 1)`.
- A window is "bad" when `b > slope_crit` (default `-0.5`) — i.e., the slope is not steep enough negative.
- Interpolate the per-window bad flag to sample rate with `interp1(..., 'nearest', ...)`. NaN tails (before first / after last window center) default to bad.

**1e. Gap interpolation** — any of the above mark a sample `bad_inds(i) = true`. Linear interpolation fills those holes into a separate `data_fixed` array so Stages 2–3 see a continuous signal without filter transients from NaN-filled gaps.

---

### Stage 2 & 3 — Band-specific envelope detection

The same iterative pipeline runs twice, once per frequency band, with different high-pass cutoffs and thresholds. Implemented in the local `compute_artifacts(...)` helper.

**Filter design** — zero-phase 4th-order Butterworth IIR high-pass, designed via `designfilt('highpassiir', ...)` with `PassbandRipple=0.2`. Two variants by default:
- **High-frequency (HF)**: cutoff `hf_pass = 35 Hz`, threshold `hf_crit = 5.5`.
  Targets muscle (EMG) contamination and electrical noise.
- **Broadband (BB)**: cutoff `bb_pass = 0.1 Hz`, threshold `bb_crit = 5.5`.
  Targets slow drifts and large-amplitude motion/sweat artifacts that survive DC removal.

Filters can be prebuilt and passed via `'hpFilt_high'` / `'hpFilt_broad'` to avoid re-designing on repeated calls (e.g., in a per-subject loop). Passing `'return_filts_only', true` returns `[hpFilt_high, hpFilt_broad]` and exits — useful for caching.

**Per-band pipeline:**

1. **Zero-phase filter** the full `data_fixed` with `filtfilt`. This removes group delay so artifact timing is accurate.
2. **Envelope**: `abs(hilbert(y))` — instantaneous amplitude.
3. **Smooth**: `movmean` over `smooth_duration = 2 s`.
4. **Log transform**: `log(y)` — converts the approximately log-normal envelope to an approximately normal distribution, so z-score statistics are meaningful.
5. **Detrend** *(on by default)*: subtract `movmedian(y, detrend_duration = 300 s)`. Moving median is chosen over mean because it's robust to the outliers we're about to detect — using the mean would let big artifacts pull the trend toward themselves. 300 s (5 min) is a good compromise: short enough to track sleep-stage-related envelope changes, long enough to be insensitive to transient spikes.
6. **Iterative z-score**:
   - Seed `detected_artifacts = bad_inds` (Stage 1 flags).
   - Compute `median` and `MAD` (or `mean`/`std` if `zscore_method='standard'`) over `~detected_artifacts & ~isexcluded`.
   - Z-score the whole envelope: `z = (y - μ) / σ`.
   - Flag `|z| > crit`.
   - **Loop**: recompute statistics from the remaining unflagged samples, re-z-score, re-threshold. Repeat until no new samples exceed the criterion.

The iterative refinement is the key idea — a single pass would have artifact energy polluting the statistics. Each iteration purges the current outliers, tightens the distribution, and exposes the next tier of artifacts. Typical recordings converge in 3–6 iterations.

**Why robust statistics (default)** — with `zscore_method='robust'`, the median is ~0.77× efficient relative to the mean for Gaussian data, but it is not dragged around by outliers. MAD (median absolute deviation) is similarly robust vs std. For sleep EEG, where 1–10% of a full night can be artifact, `robust` converges to the correct distribution faster and is the lab-standard default.

---

### Stage 4 — Combine, optionally buffer

**Union:**
```
artifacts = hf_artifacts | bb_artifacts | bad_inds
```

Any single criterion is sufficient to flag a sample.

**Exclusion mode** (`exclude_mode`):
- `'data'` (default) — samples marked `isexcluded` are ignored when computing z-score statistics, but they are *not* forced into the output. This is for cases where you want to suppress (e.g.) a known stimulus period from biasing the z-score without classifying it as artifact per se.
- `'artifact'` — the `isexcluded` mask is OR'd into the output. Use this when you want the excluded region flagged as artifact downstream.

**Buffer** (`buffer_duration`, default 0) — extend each contiguous artifact run by ±`buffer_duration` seconds on both sides. Detected via `consecutive_runs(artifacts)`, then each run's start/end indices get stretched and the buffer block is set true. Useful when filtering or spectrogram analysis can have edge effects up to a window half-length around a real artifact.

---

## Defaults at a glance

| Parameter | Default | What it controls |
|---|---|---|
| `zscore_method` | `'robust'` | median/MAD vs mean/std |
| `hf_crit` | 5.5 | HF z-score threshold |
| `hf_pass` | 35 Hz | HF high-pass cutoff |
| `hf_detrend` | `true` | remove moving median from HF envelope |
| `bb_crit` | 5.5 | BB z-score threshold |
| `bb_pass` | 0.1 Hz | BB high-pass cutoff |
| `bb_detrend` | `true` | remove moving median from BB envelope |
| `slope_test` | `true` | enable 1/f slope check |
| `slope_crit` | -0.5 | min acceptable log-log slope |
| `smooth_duration` | 2 s | envelope smoothing window |
| `detrend_duration` | 300 s | moving-median detrend window |
| `buffer_duration` | 0 s | edge buffer around detected artifacts |

### Legacy / pre-lab-standard settings

Earlier versions used `zscore_method='standard'`, `hf_crit=bb_crit=4.5`, `slope_test=false`. To reproduce results from that era:

```matlab
artifacts_old = detect_artifacts(data, Fs, ...
    'zscore_method', 'standard', ...
    'hf_crit', 4.5, 'bb_crit', 4.5, ...
    'slope_test', false);
```

---

## Diagnostic modes

- `'verbose'` → prints parameter summary and per-band iteration counts.
- `'diagnostic_plot'` → opens a stacked-axes figure per band showing filtered signal, Hilbert magnitude, smoothed magnitude, log-smoothed magnitude, and (if `*_detrend=true`) the detrended log-smoothed magnitude, all linked on x.
- `'histogram_plot'` → shows the z-score histogram redrawing at each iteration of the iterative loop. Useful to visualize convergence.

---

## Output shape

`artifacts` is always returned as a `1xN` row vector (`N = length(data)`), regardless of whether `data` was row or column. Use `~artifacts` to index the good samples.

---

## Example

```matlab
Fs = 200;                           % lab-standard sleep EEG rate
data = load_my_eeg(subject_id);     % Nx1 double

% Basic call — robust defaults, both bands, slope test on
artifacts = detect_artifacts(data, Fs);

% Verbose with a 1-second buffer around each artifact
artifacts = detect_artifacts(data, Fs, ...
    'verbose', true, ...
    'buffer_duration', 1);

% Cache filters for per-epoch looping
filts = detect_artifacts([], Fs, 'return_filts_only', true);
hpFilt_high  = filts(1);
hpFilt_broad = filts(2);
for i = 1:n_epochs
    ep_artifacts = detect_artifacts(epochs{i}, Fs, ...
        'hpFilt_high', hpFilt_high, ...
        'hpFilt_broad', hpFilt_broad);
    ...
end
```

---

## Dependencies

- `multitaper_spectrogram_mex` — from the lab's `multitaper` submodule (required when `slope_test=true`)
- `get_chunks` — from `utils/data_processing/`
- `consecutive_runs` — from `utils/data_processing/` (required when `buffer_duration > 0`)
- Signal Processing Toolbox: `designfilt`, `filtfilt`, `hilbert`, `polyfit`, `movmean`, `movmedian`

All lab dependencies are committed as siblings in this repo; `get_chunks` and `consecutive_runs` are small, self-contained helpers.

---

## Files in this directory

- `detect_artifacts.m` — the entry point and all supporting local functions
- `README.md` — this file
