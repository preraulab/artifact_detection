// Stylistic clippy lints suppressed crate-wide. Index-based loops are
// clearer than iterator chains for the numerical code (multiple parallel
// arrays, fixed-step kernels). The doc-overindentation lint trips on
// ASCII-art math indentation that mirrors the scipy reference.
#![allow(
    clippy::doc_overindented_list_items,
    clippy::manual_abs_diff,
    clippy::manual_is_multiple_of,
    clippy::manual_memcpy,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
)]

//! Rust port of `detect_artifacts.m` (Prerau Lab artifact_detection toolbox).
//!
//! Two-band envelope detection (broadband 0.1 Hz HP + high-frequency 35 Hz HP)
//! with iterative robust z-score outlier flagging, plus a sample-level
//! flat-run mask and an "outlier noise" 10-sigma backstop. Validated against
//! MATLAB R2025a `detect_artifacts(data, Fs)` on a 300 s NREM segment and on
//! a synthetic-EEG test battery (`tests/matlab_parity.rs`).
//!
//! The optional `slope_test` feature pulls in `multitaper_rs` to add the
//! MTS + log-log polyfit branch from MATLAB's `slope_test=true` path. Without
//! the feature, only the band-detection path is built.
//!
//! Public entry points:
//! - [`detect_artifacts`] — band detection on the raw signal.
//! - [`ArtifactOpts`] / [`ZScoreMethod`] — parameters.
//! - [`flat_run_mask`] / [`iterative_robust_zscore`] — exposed for tests.

pub mod bands;
pub mod filter_design;
pub mod signal;

#[cfg(feature = "slope_test")]
pub mod slope_test;

pub use bands::{
    detect_artifacts, flat_run_mask, iterative_robust_zscore, ArtifactOpts, ZScoreMethod,
};

#[cfg(feature = "slope_test")]
pub use slope_test::{detect_artifacts_with_slope, slope_test_artifacts, SlopeOpts};
