# Test-Driven Verification Record

## User journeys

- A researcher can cancel a known narrow-band component from one or more channels while retaining an unrelated component.
- A researcher can estimate a carrier and significant harmonics from multichannel data.
- A researcher can opt into spectral-amplitude matching only when explicitly requested and understand that it changes the time series.
- A researcher can inspect command-line defaults without requiring access to a recording.

## RED / GREEN evidence

- RED: `python -m pytest -q tests/test_signal.py` failed during collection because `tacs_artifact_removal` did not yet exist.
- GREEN: `python -m pytest -q` passed: `7 passed`.
- Compilation: `python -m compileall -q src tests` completed successfully.
- CLI: `PYTHONPATH=src python -m tacs_artifact_removal --help` displayed the command options, including opt-in background matching.

## Guarantees and limits

| Behavior | Test | Result |
| --- | --- | --- |
| Parabolic peak interpolation refines a sub-bin peak | `test_refine_peak_near_uses_sub_bin_interpolation` | Pass |
| Sliding least squares suppresses a target tone while preserving an unrelated tone | `test_sliding_ls_suppresses_target_tone_and_preserves_other_tone` | Pass |
| Channel-first arrays and records shorter than the requested window are handled | `test_sliding_ls_accepts_channel_first_data_and_short_records` | Pass |
| Frequencies at Nyquist are rejected | `test_sliding_ls_rejects_frequencies_at_or_above_nyquist` | Pass |
| FFT background matching changes a notched time series | `test_background_match_is_explicit_and_changes_time_series` | Pass |
| Carrier and harmonic are detected from synthetic multichannel data | `test_frequency_detection_finds_carrier_and_harmonic` | Pass |
| Background matching is disabled by default in the CLI | `test_cli_keeps_time_domain_background_matching_opt_in` | Pass |

Coverage was not measured because `pytest-cov` and `coverage` are not installed in the available Python environment. EDF read/write integration was not exercised here; the installed MNE configuration is outside the writable workspace. Validate the CLI against a local copy of an EDF recording before research use.
