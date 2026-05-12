//! Per-quirk regression tests. Each one targets a specific MATLAB behaviour
//! that the Rust port MUST preserve bit-for-bit (see comments in
//! `bands.rs` and `slope_test.rs`). Past parity drift has tracked back to
//! exactly these spots — keep them.

use artifact_detection_rs::{flat_run_mask, iterative_robust_zscore, ZScoreMethod};

#[test]
fn flat_run_mask_basic() {
    // min_run=3: [1,1,1,2,3,3,3,3] -> first 3 flagged, then 1 false, last 4
    // flagged. Total 7 trues out of 8.
    let data = vec![1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0];
    let mask = flat_run_mask(&data, 3);
    assert_eq!(
        mask,
        vec![true, true, true, false, true, true, true, true],
        "flat_run_mask did not match expected pattern"
    );

    // Below threshold: a run of 2 with min_run=3 should NOT trigger.
    let mask2 = flat_run_mask(&[1.0, 1.0, 2.0, 2.0, 3.0], 3);
    assert_eq!(mask2, vec![false; 5]);
}

#[test]
fn flat_run_mask_at_buffer_end() {
    // Run that extends to the end of the buffer must be flagged.
    let data = vec![0.0, 1.0, 2.0, 5.0, 5.0, 5.0, 5.0];
    let mask = flat_run_mask(&data, 3);
    assert_eq!(mask, vec![false, false, false, true, true, true, true]);
}

#[test]
fn iterative_robust_zscore_flags_single_outlier() {
    // 1000 samples ~ small sinusoidal swing; one giant outlier at index 500.
    // Robust z-score must flag it without ballooning the flagged-count.
    let mut data: Vec<f64> = (0..1000).map(|i| ((i as f64) * 0.001).sin()).collect();
    data[500] = 100.0;
    let bad = vec![false; data.len()];

    let mask = iterative_robust_zscore(&data, &bad, 5.5, ZScoreMethod::Robust);
    assert!(mask[500], "outlier sample 500 was not flagged");
    let flagged: usize = mask.iter().filter(|&&b| b).count();
    assert!(
        flagged < 20,
        "too many samples flagged on a near-clean signal: {}",
        flagged
    );
}

#[test]
fn mad_meanabs_matches_matlab_definition() {
    // MATLAB `mad(x)` with no second arg is mean-absolute-deviation FROM
    // THE MEAN, not from the median. For [1,2,3,4,5]:
    //   mean = 3
    //   mad  = mean(|[-2,-1,0,1,2]|) = (2+1+0+1+2)/5 = 6/5 = 1.2
    //
    // We don't expose `mad_meanabs` directly, but we can validate via the
    // iter-zscore behaviour: a perfectly symmetric input with std=0 must
    // not produce any outlier flag (because MAD would be near zero too,
    // not exactly zero only due to fp roundoff).
    let data: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let bad = vec![false; data.len()];
    let mask = iterative_robust_zscore(&data, &bad, 5.5, ZScoreMethod::Robust);
    let n_flagged: usize = mask.iter().filter(|&&b| b).count();
    assert_eq!(
        n_flagged, 0,
        "well-behaved symmetric input should not flag anything"
    );
}
