"""Numerical signal operations used by the artifact-cancellation workflow."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def _validate_signal(data: np.ndarray) -> tuple[np.ndarray, bool]:
    values = np.asarray(data, dtype=float)
    was_1d = values.ndim == 1
    if was_1d:
        values = values[np.newaxis, :]
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0:
        raise ValueError("data must be a non-empty 1D or channel-first 2D array")
    if not np.isfinite(values).all():
        raise ValueError("data must contain only finite values")
    return values, was_1d


def _validate_frequencies(
    frequencies_hz: Iterable[float], sfreq: float
) -> tuple[float, ...]:
    frequencies = tuple(float(value) for value in frequencies_hz)
    if not np.isfinite(sfreq) or sfreq <= 0:
        raise ValueError("sfreq must be a positive finite number")
    if not frequencies:
        raise ValueError("at least one artifact frequency is required")
    if any(not np.isfinite(value) or value <= 0 or value >= sfreq / 2 for value in frequencies):
        raise ValueError("artifact frequencies must be finite and strictly below Nyquist")
    return frequencies


def refine_peak_near(
    frequencies: np.ndarray,
    power: np.ndarray,
    target_hz: float,
    search_hz: float = 2.0,
) -> float:
    """Find a local spectral maximum and refine its location by parabolic interpolation."""
    frequencies = np.asarray(frequencies, dtype=float)
    power = np.asarray(power, dtype=float)
    if frequencies.ndim != 1 or power.ndim != 1 or frequencies.size != power.size:
        raise ValueError("frequencies and power must be equal-length 1D arrays")
    if frequencies.size < 3 or np.any(np.diff(frequencies) <= 0):
        raise ValueError("at least three strictly increasing frequency bins are required")
    if not np.isfinite(frequencies).all() or not np.isfinite(power).all():
        raise ValueError("frequency and power arrays must be finite")
    if not np.isfinite(target_hz) or not np.isfinite(search_hz) or search_hz <= 0:
        raise ValueError("target_hz must be finite and search_hz must be positive")

    indices = np.flatnonzero(np.abs(frequencies - target_hz) <= search_hz)
    if indices.size == 0:
        return float(target_hz)
    index = int(indices[np.argmax(power[indices])])
    if index == 0 or index == power.size - 1:
        return float(frequencies[index])

    left, center, right = power[index - 1 : index + 2]
    denominator = left - 2.0 * center + right
    if denominator == 0:
        return float(frequencies[index])
    offset = 0.5 * (left - right) / denominator
    bin_width = frequencies[index + 1] - frequencies[index]
    return float(frequencies[index] + offset * bin_width)


def _build_reference(
    n_samples: int, sfreq: float, frequencies_hz: tuple[float, ...]
) -> np.ndarray:
    time = np.arange(n_samples, dtype=float) / sfreq
    columns: list[np.ndarray] = []
    for frequency in frequencies_hz:
        phase = 2.0 * np.pi * frequency * time
        columns.extend((np.sin(phase), np.cos(phase)))
    return np.column_stack(columns)


def sliding_ls_cancel(
    data: np.ndarray,
    sfreq: float,
    frequencies_hz: Iterable[float],
    window_seconds: float = 1.0,
    hop_seconds: float = 0.1,
) -> np.ndarray:
    """Subtract sinusoidal components estimated independently in overlapping windows.

    Accepts a single time series or a channel-first ``(channels, samples)`` array.
    The shape of the input is preserved. Frequency selection is shared across
    channels, while each channel gets its own least-squares coefficients.
    """
    values, was_1d = _validate_signal(data)
    frequencies = _validate_frequencies(frequencies_hz, sfreq)
    if not np.isfinite(window_seconds) or window_seconds <= 0:
        raise ValueError("window_seconds must be a positive finite number")
    if not np.isfinite(hop_seconds) or hop_seconds <= 0:
        raise ValueError("hop_seconds must be a positive finite number")

    n_samples = values.shape[1]
    window = min(n_samples, max(1, int(round(window_seconds * sfreq))))
    hop = max(1, int(round(hop_seconds * sfreq)))
    n_basis = 2 * len(frequencies)
    if window < n_basis:
        raise ValueError("each fitting window must contain at least two samples per frequency")

    reference = _build_reference(n_samples, sfreq, frequencies)
    starts = list(range(0, n_samples - window + 1, hop))
    final_start = n_samples - window
    if not starts or starts[-1] != final_start:
        starts.append(final_start)

    cleaned = np.empty_like(values)
    for channel_index, channel in enumerate(values):
        estimate_sum = np.zeros(n_samples, dtype=float)
        coverage = np.zeros(n_samples, dtype=float)
        for start in starts:
            stop = start + window
            window_reference = reference[start:stop]
            coefficients, *_ = np.linalg.lstsq(
                window_reference, channel[start:stop], rcond=None
            )
            estimate_sum[start:stop] += window_reference @ coefficients
            coverage[start:stop] += 1.0
        cleaned[channel_index] = channel - estimate_sum / coverage

    return cleaned[0] if was_1d else cleaned


def background_match_fft(
    data: np.ndarray,
    sfreq: float,
    frequencies_hz: Iterable[float],
    half_width_hz: float = 1.5,
    background_width_hz: float = 3.0,
) -> np.ndarray:
    """Raise selected FFT-bin amplitudes to a local background floor.

    This operation is intentionally opt-in: inverse-transforming the modified
    spectrum changes the time series and is not a lossless visualization step.
    """
    values = np.asarray(data, dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.isfinite(values).all():
        raise ValueError("data must be a non-empty finite 1D array")
    frequencies = _validate_frequencies(frequencies_hz, sfreq)
    if not np.isfinite(half_width_hz) or half_width_hz <= 0:
        raise ValueError("half_width_hz must be positive and finite")
    if not np.isfinite(background_width_hz) or background_width_hz <= 0:
        raise ValueError("background_width_hz must be positive and finite")

    spectrum = np.fft.rfft(values)
    bins = np.fft.rfftfreq(values.size, d=1.0 / sfreq)
    for frequency in frequencies:
        distance = np.abs(bins - frequency)
        target = distance <= half_width_hz
        background = (distance > half_width_hz) & (
            distance <= half_width_hz + background_width_hz
        )
        if np.count_nonzero(background) < 10:
            continue
        floor = np.median(np.abs(spectrum[background]))
        target_indices = np.flatnonzero(target)
        below_floor = np.abs(spectrum[target_indices]) < floor
        selected = target_indices[below_floor]
        spectrum[selected] = floor * np.exp(1j * np.angle(spectrum[selected]))

    return np.fft.irfft(spectrum, n=values.size)
