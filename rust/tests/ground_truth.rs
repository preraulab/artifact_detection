//! Ground-truth recovery tests — MATLAB-independent.
//!
//! These tests assert that the Rust artifact detector
//! 1. flags at least one sample inside a tight window around each injected
//!    motion-artifact center (per-event recovery), and
//! 2. does not flag a non-trivial fraction of a fully-clean signal.
//!
//! They run on the synthetic-EEG fixtures from
//! `DYNAM-O_DesktopApp/scripts/generate_handoff_fixtures.m` — no MATLAB
//! reference output is loaded, so this suite catches regressions even on
//! CI runners without MATLAB.

use artifact_detection_rs::{detect_artifacts, ArtifactOpts};
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

/// Parse the motion-event center times (seconds) from the JSON sidecar
/// emitted by `generate_handoff_fixtures.m`.
fn motion_times_from_events(variant: &str, seed: u32) -> Vec<f64> {
    let path = fixtures().join(format!("sim_{}_seed{}_events.json", variant, seed));
    let json = std::fs::read_to_string(&path).expect("read events.json");
    // jsonencode -> '"motion_times_s": [t1, t2, ...]' (might be on one line)
    let key = "\"motion_times_s\":";
    let i = json.find(key).expect("motion_times_s key missing");
    let after = &json[i + key.len()..];
    let lb = after.find('[').expect("missing [");
    let rb = after.find(']').expect("missing ]");
    after[lb + 1..rb]
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect()
}

#[test]
fn sim_motion_recovers_each_event() {
    let data = load_double(fixtures().join("sim_motion_seed42.npy"));
    let fs = 100.0_f64;
    let opts = ArtifactOpts::default();
    let mask = detect_artifacts(&data, fs, &opts);

    let centers = motion_times_from_events("motion", 42);
    assert!(!centers.is_empty(), "no motion centers found in events.json");

    let tolerance_samples: usize = (0.2 * fs) as usize; // +/-200 ms recovery window
    for c in &centers {
        let center_idx = (*c * fs).round() as i64;
        let lo = (center_idx - tolerance_samples as i64).max(0) as usize;
        let hi = ((center_idx + tolerance_samples as i64) as usize).min(mask.len());
        let hit = mask[lo..hi].iter().any(|&b| b);
        assert!(
            hit,
            "motion event at t={} s (idx={}) not recovered within +/-{} samples",
            c, center_idx, tolerance_samples
        );
    }
}

#[test]
fn sim_clean_low_false_positive_rate() {
    let data = load_double(fixtures().join("sim_clean_seed42.npy"));
    let fs = 100.0_f64;
    let opts = ArtifactOpts::default();
    let mask = detect_artifacts(&data, fs, &opts);

    let n = mask.len();
    let flagged: usize = mask.iter().filter(|&&b| b).count();
    let pct = flagged as f64 / n as f64;
    // Band-only detection (no slope test edge zone) on a clean record should
    // be well under 1% false-positive.
    assert!(
        pct < 0.01,
        "clean signal flagged {}/{} = {:.3}% (>1%)",
        flagged, n, 100.0 * pct
    );
}

#[test]
fn sim_flat_covers_injected_runs() {
    let data = load_double(fixtures().join("sim_flat_seed42.npy"));
    let fs = 100.0_f64;
    let opts = ArtifactOpts::default();
    let mask = detect_artifacts(&data, fs, &opts);

    // Injected flat intervals (1-based, inclusive from MATLAB):
    //   [6001 .. 6151], [22501 .. 22651]
    // Convert to 0-based half-open.
    for (a_mat, b_mat) in [(6001usize, 6151usize), (22501, 22651)] {
        let a = a_mat - 1;
        let b = b_mat - 1;
        let covered = mask[a..=b].iter().filter(|&&x| x).count();
        let len = b - a + 1;
        let frac = covered as f64 / len as f64;
        assert!(
            frac > 0.95,
            "flat run [{}..{}] only {:.1}% covered (<95%)",
            a, b, 100.0 * frac
        );
    }
}
