<div align="center">
<br>

<h1>tACS / Hi-tACS Artifact Cancellation</h1>

<p><strong>A frequency-informed signal-processing workflow for concurrent EEG / SEEG.</strong></p>
<p>Sliding-window sinusoidal least squares with explicit assumptions and analysis boundaries.</p>
</div>

![Representative single-channel EEG comparison before and after processing: time series, power spectral density, and short-time Fourier transform](assets/result-comparison.png)

<p align="center"><sub>Illustrative source result only. The plotted output includes least-squares subtraction and the original script's spectral-amplitude modification; it is not an isolated validation of either operation.</sub></p>

---

## Computational Model

The workflow estimates a shared set of candidate frequencies from the channel-averaged Welch power spectral density (PSD). For each channel, sine/cosine reference coefficients are fitted independently within overlapping windows. The fitted component is subtracted from the observation.

This is a frequency-informed linear projection, not a classifier that distinguishes stimulation artifact from neural activity. Any neural component that overlaps the selected reference subspace can also be attenuated. A reduction in power at a target frequency is therefore not, by itself, evidence of artifact-specific removal or neural-signal preservation.

## Processing Stages

| Stage | Operation | Main assumption |
| --- | --- | --- |
| Frequency estimation | Refine the nominal carrier; screen harmonics and optional sidebands using local PSD ratios | Artifact frequencies are known approximately and produce detectable narrow-band peaks |
| Local fit | Estimate channel-specific sine/cosine coefficients by least squares in overlapping windows | Frequency is approximately fixed within each window; amplitude and phase may vary between windows |
| Subtraction | Average overlapping fitted waveforms and subtract from the observation | The modeled component is a useful approximation of the artifact contribution |
| Assessment | Inspect residual time series, PSD, and STFT; evaluate signal retention separately | Artifact suppression and neural-signal preservation require distinct evidence |

## Installation and Use

Python 3.10 or newer.

    pip install -e ".[edf,plot]"
    tacs-artifact-removal recording.edf --frequency 40 --plot

The command reads an EDF file, estimates candidate frequencies, performs sliding-window least-squares subtraction, and writes recording_cleaned_raw.fif. The nominal carrier is required input. Use --help for all options.

Sideband screening is opt-in with --include-sidebands and should be reviewed against the stimulation protocol and channel-level spectra.

## Spectral-Amplitude Modification

The optional --background-match operation raises selected FFT-bin amplitudes to a local background floor and then applies an inverse FFT. It changes the time-domain signal as well as its spectrum; preserving the phase angle of modified bins does not preserve the original waveform. The operation is disabled by default and should not be used for power, phase, connectivity, or other quantitative analyses unless its consequences are the subject of a separately validated analysis.

## Scientific Scope

The model targets additive narrow-band components with approximately stable frequencies and locally varying amplitudes/phases. It does not explicitly track continuous frequency drift, model stimulation onset/offset transients, recover clipped data, or establish whether a fitted component is neural or artifactual. Neural activity at the modeled frequencies may be removed.

The included figure is an example from the original analysis script, which applied both least-squares subtraction and spectral-amplitude modification before plotting. It cannot identify the contribution of either step or establish neural-signal fidelity, generalization across participants, or clinical validity. The method document details assumptions, numerical behavior, and a validation framework.

## Repository

| Path | Contents |
| --- | --- |
| src/tacs_artifact_removal/ | Reusable signal processing, frequency estimation, and EDF command line |
| docs/methodology.md | Signal model, numerical method, assumptions, and validation considerations |
| tests/ | Synthetic-signal regression tests |
| assets/result-comparison.png | Original example figure |

Install the development extra and run the tests with:

    pip install -e ".[dev]"
    pytest

## License

No license is included. Reuse and redistribution permissions should not be assumed until the author selects a license.
