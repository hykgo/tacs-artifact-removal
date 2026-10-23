"""Sliding-window sinusoidal artifact cancellation for EEG/SEEG."""

from .frequencies import DetectionResult, detect_artifact_frequencies
from .signal import background_match_fft, refine_peak_near, sliding_ls_cancel

__all__ = [
    "DetectionResult",
    "background_match_fft",
    "detect_artifact_frequencies",
    "refine_peak_near",
    "sliding_ls_cancel",
]

__version__ = "0.1.0"
