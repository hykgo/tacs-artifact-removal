"""Command-line interface for EDF input and FIF output."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

from .frequencies import detect_artifact_frequencies
from .signal import background_match_fft, sliding_ls_cancel


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tacs-artifact-removal",
        description="Estimate and subtract narrow-band tACS artifacts from an EDF recording.",
    )
    parser.add_argument("input", type=Path, help="input EDF file")
    parser.add_argument("--frequency", type=float, required=True, help="nominal carrier frequency in Hz")
    parser.add_argument("--output", type=Path, help="output FIF path ending in _raw.fif (default: <input>_cleaned_raw.fif)")
    parser.add_argument("--window", type=float, default=1.0, help="least-squares window in seconds (default: 1.0)")
    parser.add_argument("--hop", type=float, default=0.1, help="window hop in seconds (default: 0.1)")
    parser.add_argument("--harmonics", type=int, default=7, help="maximum harmonic to inspect (default: 7)")
    parser.add_argument("--ratio", type=float, default=10.0, help="minimum peak/background ratio (default: 10)")
    parser.add_argument("--include-sidebands", action="store_true", help="also inspect significant nearby spectral peaks")
    parser.add_argument("--background-match", action="store_true", help="apply an additional FFT amplitude-floor operation that changes the time series")
    parser.add_argument("--channel", type=int, default=0, help="channel index for optional diagnostic plot (default: 0)")
    parser.add_argument("--plot", action="store_true", help="show before/after time series and PSD for one channel")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.input.is_file():
        parser.error(f"input file does not exist: {args.input}")
    if args.output is not None and not args.output.name.lower().endswith("_raw.fif"):
        parser.error("--output must end in _raw.fif")

    try:
        import mne
    except ImportError as error:
        print("EDF I/O requires the optional dependency: pip install .[edf]", file=sys.stderr)
        raise SystemExit(2) from error

    raw = mne.io.read_raw_edf(args.input, preload=True, verbose="error")
    sfreq = float(raw.info["sfreq"])
    data = raw.get_data()
    detection = detect_artifact_frequencies(
        data,
        sfreq=sfreq,
        nominal_frequency=args.frequency,
        n_harmonics=args.harmonics,
        peak_background_ratio=args.ratio,
        include_sidebands=args.include_sidebands,
    )
    if not detection.frequencies_hz:
        parser.error("no significant artifact frequencies detected; inspect the PSD and parameters")

    cleaned = sliding_ls_cancel(
        data,
        sfreq=sfreq,
        frequencies_hz=detection.frequencies_hz,
        window_seconds=args.window,
        hop_seconds=args.hop,
    )
    if args.background_match:
        print("Warning: background matching modifies the time-domain signal.", file=sys.stderr)
        cleaned = np.stack(
            [background_match_fft(channel, sfreq, detection.frequencies_hz) for channel in cleaned]
        )

    output = args.output or args.input.with_name(f"{args.input.stem}_cleaned_raw.fif")
    cleaned_raw = mne.io.RawArray(
        cleaned,
        raw.info.copy(),
        first_samp=raw.first_samp,
        verbose="error",
    )
    cleaned_raw.set_annotations(raw.annotations.copy())
    cleaned_raw.save(output, overwrite=True, verbose="error")

    print(f"Sampling rate: {sfreq:g} Hz; channels: {len(raw.ch_names)}")
    print(f"Carrier estimate: {detection.carrier_hz:.4f} Hz")
    print("Cancellation frequencies: " + ", ".join(f"{value:.4f}" for value in detection.frequencies_hz))
    print(f"Saved: {output}")

    if args.plot:
        if not 0 <= args.channel < data.shape[0]:
            parser.error(f"--channel must be between 0 and {data.shape[0] - 1}")
        _show_diagnostic(data[args.channel], cleaned[args.channel], sfreq, args.channel)
    return 0


def _show_diagnostic(before, after, sfreq: float, channel: int) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.signal import welch

    time = np.arange(before.size) / sfreq
    nperseg = min(before.size, max(8, int(sfreq * 8)))
    frequencies, before_psd = welch(before, fs=sfreq, nperseg=nperseg)
    _, after_psd = welch(after, fs=sfreq, nperseg=nperseg)
    figure, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    axes[0].plot(time, before, linewidth=0.6, label="Before")
    axes[0].plot(time, after, linewidth=0.6, label="After")
    axes[0].set(title=f"Channel {channel}: time domain", xlabel="Time (s)", ylabel="Amplitude")
    axes[0].legend(frameon=False)
    axes[1].semilogy(frequencies, before_psd, label="Before")
    axes[1].semilogy(frequencies, after_psd, label="After")
    axes[1].set(title="Welch PSD", xlabel="Frequency (Hz)", ylabel="Power / Hz")
    axes[1].legend(frameon=False)
    plt.show()
