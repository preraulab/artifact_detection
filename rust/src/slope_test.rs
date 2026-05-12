//! Slope-test artifact detection — port of `detect_artifacts.m::slope_test`
//! (lines 142–161). Computes a multitaper spectrogram with `[NW=10, K=19]`
//! on 10 s / 5 s-step windows over `[1, min(55, Fs/2)] Hz`, fits a linear
//! regression of `log(power) vs log(freq)` per window, and flags samples
//! whose window slope exceeds `slope_crit` (default `-0.5`).
//!
//! # MATLAB quirk preserved
//!
//! Per `detect_artifacts.m` lines 156–158:
//! ```matlab
//! bad_slope = interp1(stimes, double(B(:,1)>slope_crit)', t, 'nearest')';
//! bad_slope(isnan(bad_slope)) = 1;
//! ```
//! Out-of-range samples (the first and last `step_s` of the recording for a
//! 10 s / 5 s window) are flagged unconditionally — *even if* the nearest
//! window's slope verdict was "good". The Rust port keeps the same
//! semantics; missing this branch caused a ~800-sample edge under-flag bug
//! during initial porting.

use multitaper_rs::{
    compute_spectrogram, dpss, DetrendMode, DpssError, SpectrogramParams, Weighting,
};
use ndarray::Array1;

use crate::bands::{detect_artifacts as rs_detect_artifacts, ArtifactOpts, ZScoreMethod};

#[derive(Debug, thiserror::Error)]
pub enum SlopeError {
    #[error("DPSS error: {0}")]
    Dpss(#[from] DpssError),
    #[error("spectrogram error: {0}")]
    Spectrogram(String),
}

/// Options for the slope-test branch.
#[derive(Debug, Clone)]
pub struct SlopeOpts {
    pub window_s: f64,
    pub step_s: f64,
    pub nw: f64,
    pub k: usize,
    pub slope_crit: f64,
    pub fmin: f64,
    pub fmax_cap: f64,
}

impl Default for SlopeOpts {
    fn default() -> Self {
        Self {
            window_s: 10.0,
            step_s: 5.0,
            nw: 10.0,
            k: 19,
            slope_crit: -0.5,
            fmin: 1.0,
            fmax_cap: 55.0,
        }
    }
}

fn next_pow2(x: f64) -> usize {
    if x <= 1.0 {
        return 1;
    }
    let mut p = 1usize;
    while (p as f64) < x {
        p <<= 1;
    }
    p
}

/// Compute the slope-test artifact mask. Returns a length-`data.len()`
/// boolean mask. On DPSS / spectrogram failure returns `Err` and the
/// caller decides whether to treat that as all-false or propagate.
pub fn slope_test_artifacts(
    data: &[f64],
    fs: f64,
    opts: &SlopeOpts,
) -> Result<Vec<bool>, SlopeError> {
    let win_samples = (opts.window_s * fs).round() as usize;
    let fmax = opts.fmax_cap.min(fs / 2.0);
    let freq_range = (opts.fmin, fmax);
    let dsfreqs = 0.1_f64;
    let nfft = next_pow2(fs / dsfreqs);
    if win_samples < 4 || data.len() < win_samples {
        return Ok(vec![false; data.len()]);
    }

    let (tapers, _ratios) = dpss(win_samples, opts.nw, opts.k)?;
    let out = compute_spectrogram(
        Array1::from(data.to_vec()).view(),
        tapers.view(),
        None,
        &SpectrogramParams {
            fs,
            frequency_range: freq_range,
            window_params: (opts.window_s, opts.step_s),
            nfft,
            detrend: DetrendMode::Constant,
            weighting: Weighting::Unity,
        },
    )
    .map_err(SlopeError::Spectrogram)?;

    let spect = &out.mt_spectrogram; // (F, T_win)
    let stimes = &out.stimes;
    let sfreqs = &out.sfreqs;
    let nf = sfreqs.len();
    let nt = stimes.len();
    if nf == 0 || nt == 0 {
        return Ok(vec![false; data.len()]);
    }

    let log_f: Vec<f64> = sfreqs.iter().map(|f| f.max(1e-30).ln()).collect();
    let mean_log_f: f64 = log_f.iter().sum::<f64>() / nf as f64;
    let denom: f64 = log_f.iter().map(|x| (x - mean_log_f).powi(2)).sum();
    if denom <= 0.0 || !denom.is_finite() {
        return Ok(vec![false; data.len()]);
    }

    // Per-window slope sign test.
    let mut bad_win: Vec<bool> = vec![false; nt];
    for t in 0..nt {
        let mut log_s: Vec<f64> = Vec::with_capacity(nf);
        let mut any_bad = false;
        for f in 0..nf {
            let v = spect[[f, t]];
            if !(v > 0.0 && v.is_finite()) {
                any_bad = true;
                break;
            }
            log_s.push(v.ln());
        }
        if any_bad {
            bad_win[t] = true;
            continue;
        }
        let mean_log_s: f64 = log_s.iter().sum::<f64>() / nf as f64;
        let mut num = 0.0_f64;
        for f in 0..nf {
            num += (log_f[f] - mean_log_f) * (log_s[f] - mean_log_s);
        }
        let slope = num / denom;
        bad_win[t] = slope > opts.slope_crit;
    }

    // Per-sample nearest-window assignment with out-of-range flagged true.
    // This out-of-range branch is the documented MATLAB quirk — DO NOT remove.
    let n = data.len();
    let t_first = stimes[0];
    let t_last = stimes[nt - 1];
    let mut out_mask: Vec<bool> = vec![false; n];
    for i in 0..n {
        let t = i as f64 / fs;
        if !t.is_finite() || t < t_first || t > t_last {
            out_mask[i] = true;
            continue;
        }
        let mut lo = 0_usize;
        let mut hi = nt - 1;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if stimes[mid] < t {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let mut idx = lo;
        if idx > 0 && idx < nt {
            let d_prev = (t - stimes[idx - 1]).abs();
            let d_curr = (t - stimes[idx]).abs();
            if d_prev < d_curr {
                idx -= 1;
            }
        }
        out_mask[i] = bad_win[idx];
    }
    Ok(out_mask)
}

/// Convenience wrapper: band detection PLUS slope-test, folded in **before**
/// the band-detection iter-zscore (matches MATLAB order). The desktop-app
/// orchestrator OR's the slope mask AFTER band detection — that's a known
/// minor divergence; this function uses the MATLAB order.
pub fn detect_artifacts_with_slope(
    data: &[f64],
    fs: f64,
    band_opts: &ArtifactOpts,
    slope_opts: &SlopeOpts,
) -> Result<Vec<bool>, SlopeError> {
    let slope_mask = slope_test_artifacts(data, fs, slope_opts)?;

    // MATLAB folds bad_slope into bad_inds BEFORE the band detector runs (so
    // those samples are excluded from the iter-zscore centering stats). The
    // band detector here doesn't take an explicit `bad_inds` parameter; pre-
    // fill the slope-flagged samples with the nanmean of the rest, which is
    // the same as MATLAB's `interp1` of the unflagged signal, then run the
    // band detector and OR the slope mask back in.
    let n = data.len();
    let mut data_pre = data.to_vec();
    let good_mean = {
        let mut s = 0.0_f64;
        let mut k = 0_usize;
        for i in 0..n {
            if !slope_mask[i] {
                s += data[i];
                k += 1;
            }
        }
        if k > 0 {
            s / k as f64
        } else {
            0.0
        }
    };
    for i in 0..n {
        if slope_mask[i] {
            data_pre[i] = good_mean;
        }
    }

    let mut band_mask = rs_detect_artifacts(&data_pre, fs, band_opts);
    for i in 0..n {
        band_mask[i] = band_mask[i] || slope_mask[i];
    }
    Ok(band_mask)
}

#[doc(hidden)]
// Silence the unused-import warning if `slope_test` is built without
// downstream users of ZScoreMethod via this wrapper.
#[allow(dead_code)]
fn _hold_zscore_method(_: ZScoreMethod) {}
