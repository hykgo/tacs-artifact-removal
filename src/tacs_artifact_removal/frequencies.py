"""PSD-based detection of carrier, harmonic, and optional sideband frequencies."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import find_peaks, welch

from .signal import refine_peak_near


@dataclass(frozen=True)
class DetectionResult:
    carrier_hz: float
    frequencies_hz: tuple[float, ...]


def detect_artifact_frequencies(
    data: np.ndarray,
    sfreq: float,
    nominal_frequency: float,
    n_harmonics: int = 7,
    peak_background_ratio: float = 10.0,
    search_hz: float = 2.0,
    nperseg_seconds: float = 8.0,
    include_sidebands: bool = False,
) -> DetectionResult:
    """Estimate significant artifact frequencies from mean multichannel Welch PSD."""
    values = np.asarray(data, dtype=float)
    if values.ndim == 1:
        values = values[np.newaxis, :]
    if values.ndim != 2 or min(values.shape) == 0 or not np.isfinite(values).all():
        raise ValueError("data must be a non-empty finite 1D or channel-first 2D array")
    if not np.isfinite(sfreq) or sfreq <= 0:
        raise ValueError("sfreq must be positive and finite")
    if not np.isfinite(nominal_frequency) or not 0 < nominal_frequency < sfreq / 2:
        raise ValueError("nominal_frequency must be positive and below Nyquist")
    if not isinstance(n_harmonics, int) or n_harmonics < 1:
        raise ValueError("n_harmonics must be a positive integer")
    if not np.isfinite(peak_background_ratio) or peak_background_ratio <= 1:
        raise ValueError("peak_background_ratio must be greater than one")
    if not np.isfinite(nperseg_seconds) or nperseg_seconds <= 0:
        raise ValueError("nperseg_seconds must be positive and finite")

    nperseg = min(values.shape[1], max(3, int(round(sfreq * nperseg_seconds))))
    noverlap = min(nperseg - 1, int(round(nperseg * 0.5)))
    frequencies, channel_psd = welch(
        values, fs=sfreq, nperseg=nperseg, noverlap=noverlap, axis=-1
    )
    mean_psd = np.mean(channel_psd, axis=0)
    carrier = refine_peak_near(frequencies, mean_psd, nominal_frequency, search_hz)
    resolution = frequencies[1] - frequencies[0]
    selected: list[float] = []

    def significant_peak(target_hz: float, background_width_hz: float) -> float | None:
        if target_hz >= sfreq / 2:
            return None
        refined = refine_peak_near(frequencies, mean_psd, target_hz, search_hz)
        peak_index = int(np.argmin(np.abs(frequencies - refined)))
        distance = np.abs(frequencies - refined)
        background = mean_psd[
            (distance > max(2.0 * resolution, search_hz))
            & (distance <= background_width_hz)
        ]
        if background.size == 0:
            return None
        baseline = float(np.median(background))
        if baseline <= 0 or mean_psd[peak_index] <= baseline * peak_background_ratio:
            return None
        return float(refined)

    for harmonic in range(1, n_harmonics + 1):
        frequency = significant_peak(harmonic * carrier, background_width_hz=5.0)
        if frequency is not None:
            selected.append(frequency)

    if include_sidebands:
        peak_distance = max(1, int(round(1.0 / resolution)))
        peak_indices, _ = find_peaks(
            mean_psd,
            height=np.median(mean_psd) * peak_background_ratio,
            distance=peak_distance,
        )
        for peak_index in peak_indices:
            candidate = float(frequencies[peak_index])
            if any(abs(candidate - harmonic) <= search_hz for harmonic in selected):
                continue
            if any(2.0 < abs(candidate - harmonic) < 30.0 for harmonic in selected):
                frequency = significant_peak(candidate, background_width_hz=3.0)
                if frequency is not None:
                    selected.append(frequency)

    unique = tuple(sorted({round(frequency, 4) for frequency in selected}))
    return DetectionResult(carrier_hz=float(carrier), frequencies_hz=unique)
