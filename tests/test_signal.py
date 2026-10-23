import numpy as np
import pytest

from tacs_artifact_removal.signal import (
    background_match_fft,
    refine_peak_near,
    sliding_ls_cancel,
)
from tacs_artifact_removal.frequencies import detect_artifact_frequencies


def test_refine_peak_near_uses_sub_bin_interpolation():
    frequencies = np.arange(0.0, 21.0)
    true_peak = 10.25
    power = 1.0 - (frequencies - true_peak) ** 2

    estimate = refine_peak_near(frequencies, power, target_hz=10.0, search_hz=2.0)

    assert estimate == pytest.approx(true_peak)


def test_sliding_ls_suppresses_target_tone_and_preserves_other_tone():
    sfreq = 500.0
    time = np.arange(int(5 * sfreq)) / sfreq
    wanted = 0.4 * np.sin(2 * np.pi * 10 * time)
    artifact = 3.0 * np.sin(2 * np.pi * 40 * time + 0.4)
    observed = wanted + artifact

    cleaned = sliding_ls_cancel(
        observed,
        sfreq=sfreq,
        frequencies_hz=[40.0],
        window_seconds=1.0,
        hop_seconds=0.25,
    )

    assert np.sqrt(np.mean(cleaned**2)) < 0.5
    assert np.corrcoef(cleaned, wanted)[0, 1] > 0.98


def test_sliding_ls_accepts_channel_first_data_and_short_records():
    sfreq = 200.0
    time = np.arange(100) / sfreq
    artifact = np.vstack(
        [np.sin(2 * np.pi * 40 * time), 2 * np.cos(2 * np.pi * 40 * time)]
    )

    cleaned = sliding_ls_cancel(
        artifact,
        sfreq=sfreq,
        frequencies_hz=[40.0],
        window_seconds=1.0,
        hop_seconds=0.2,
    )

    assert cleaned.shape == artifact.shape
    assert np.max(np.abs(cleaned)) < 1e-10


def test_sliding_ls_rejects_frequencies_at_or_above_nyquist():
    with pytest.raises(ValueError, match="Nyquist"):
        sliding_ls_cancel(np.ones(100), 100.0, [50.0])


def test_background_match_is_explicit_and_changes_time_series():
    sfreq = 256.0
    rng = np.random.default_rng(4)
    spectrum = np.fft.rfft(rng.normal(size=4096))
    bins = np.fft.rfftfreq(4096, d=1 / sfreq)
    spectrum[np.abs(bins - 10.0) <= 0.5] = 0
    signal = np.fft.irfft(spectrum, n=4096)

    matched = background_match_fft(signal, sfreq, [10.0], half_width_hz=0.5)

    assert matched.shape == signal.shape
    assert not np.allclose(matched, signal)


def test_frequency_detection_finds_carrier_and_harmonic():
    sfreq = 256.0
    time = np.arange(int(12 * sfreq)) / sfreq
    data = np.vstack(
        [
            0.2 * np.sin(2 * np.pi * 12 * time)
            + 5.0 * np.sin(2 * np.pi * 40.25 * time)
            + 2.0 * np.sin(2 * np.pi * 80.5 * time),
            0.15 * np.cos(2 * np.pi * 15 * time)
            + 3.0 * np.sin(2 * np.pi * 40.25 * time + 0.3)
            + 1.0 * np.sin(2 * np.pi * 80.5 * time),
        ]
    )

    result = detect_artifact_frequencies(
        data,
        sfreq=sfreq,
        nominal_frequency=40.0,
        n_harmonics=2,
        nperseg_seconds=8.0,
    )

    assert result.carrier_hz == pytest.approx(40.25, abs=0.1)
    assert any(abs(frequency - 40.25) < 0.1 for frequency in result.frequencies_hz)
    assert any(abs(frequency - 80.5) < 0.2 for frequency in result.frequencies_hz)
