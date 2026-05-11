//! Slope-test feature tests — only compiled with `--features slope_test`.
//!
//! Verifies:
//!   1. The MATLAB edge-flag quirk is preserved: with the default
//!      window=10s / step=5s, the first and last 5 s of every record are
//!      flagged regardless of slope verdict (matches MATLAB's
//!      `bad_slope(isnan(...))=1` line). Past parity drift was ~800
//!      samples here; the test pins it at >= 95% of the edge zone.
//!   2. `detect_artifacts_with_slope` runs to completion on the
//!      synthetic-EEG sim_full fixture and produces a mask of the right
//!      length.

#![cfg(feature = "slope_test")]

use artifact_detection_rs::{
    detect_artifacts_with_slope, slope_test_artifacts, ArtifactOpts, SlopeOpts,
};
use ndarray::Array1;
use ndarray_npy::read_npy;
use std::path::PathBuf;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn load_double(path: PathBuf) -> Vec<f64> {
    let arr: Array1<f64> = read_npy(path).expect("load f64 npy");
    arr.to_vec()
}

#[test]
fn slope_test_edge_flags_5s() {
    let data = load_double(fixtures().join("sim_clean_seed42.npy"));
    let fs = 100.0_f64;
    let opts = SlopeOpts::default();
    let mask = slope_test_artifacts(&data, fs, &opts).expect("slope_test");

    let n = mask.len();
    let edge_samples: usize = (opts.step_s * fs) as usize; // 500

    // First step_s seconds: must be >= 95% flagged.
    let head: usize = mask[..edge_samples].iter().filter(|&&x| x).count();
    let head_pct = head as f64 / edge_samples as f64;
    assert!(
        head_pct >= 0.95,
        "head edge zone only {:.1}% flagged (<95%)",
        100.0 * head_pct
    );

    // Last step_s seconds: same.
    let tail: usize = mask[n - edge_samples..].iter().filter(|&&x| x).count();
    let tail_pct = tail as f64 / edge_samples as f64;
    assert!(
        tail_pct >= 0.95,
        "tail edge zone only {:.1}% flagged (<95%)",
        100.0 * tail_pct
    );
}

#[test]
fn detect_artifacts_with_slope_runs_on_full_sim() {
    let data = load_double(fixtures().join("sim_full_seed42.npy"));
    let fs = 100.0_f64;
    let band_opts = ArtifactOpts::default();
    let slope_opts = SlopeOpts::default();
    let mask = detect_artifacts_with_slope(&data, fs, &band_opts, &slope_opts)
        .expect("detect_artifacts_with_slope");
    assert_eq!(mask.len(), data.len());

    // Combined output should cover at least the injected flat-run extents.
    // (Stronger checks live in ground_truth.rs against the band-only path.)
    let flagged: usize = mask.iter().filter(|&&x| x).count();
    assert!(
        flagged > 300, // >= the 302 flat-run samples
        "combined slope+band mask flagged only {} samples on sim_full",
        flagged
    );
}
