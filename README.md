# Digital Tube Amplifier Model

A MATLAB project for data-driven digital modelling of a nonlinear tube guitar amplifier.

The project was developed as part of a bachelor’s thesis and implements a level-dependent Parallel Hammerstein model derived from measured excitation and response signals.

The model is designed to reproduce both the linear frequency response and nonlinear harmonic behaviour of a real guitar amplifier.

## Overview

A tube guitar amplifier is a nonlinear and level-dependent system. Its response changes not only with frequency, but also with the amplitude of the input signal.

To represent this behaviour, the project uses a Parallel Hammerstein structure consisting of multiple nonlinear branches:

```text
Input signal
     │
     ├── x¹[n] ── h₁[n] ──┐
     ├── x²[n] ── h₂[n] ──┤
     ├── x³[n] ── h₃[n] ──┤
     │         ...         ├── Sum ── Output signal
     └── xᴹ[n] ── hᴹ[n] ──┘
```

Each branch applies a polynomial nonlinearity to the input signal and then filters the result using an identified Hammerstein kernel.

The output is calculated as:

```text
y[n] = Σₘ hₘ[n] * xᵐ[n]
```

where:

* `x[n]` is the input signal,
* `m` is the nonlinear order,
* `hₘ[n]` is the impulse response associated with the corresponding nonlinear branch,
* `*` denotes convolution.

## Main Features

* nonlinear guitar amplifier modelling in MATLAB
* Parallel Hammerstein model architecture
* level-dependent kernel interpolation
* static and time-varying processing modes
* oversampling of nonlinear branches
* anti-aliasing during downsampling
* gain calibration against measured amplifier responses
* comparison between measured and modelled output
* WAV input and output support
* support for measured exponential sine sweep data

## Level-Dependent Model

The amplifier was measured at several input levels in dBFS.

For each measured level, a set of Parallel Hammerstein kernels is stored in the model. During processing, the program determines the current signal level and obtains a corresponding kernel set.

When the requested input level lies between two measured levels, the kernels are interpolated along the level axis.

```text
Measured kernel set at level A
               │
               ├── interpolation ── kernel set for current level
               │
Measured kernel set at level B
```

This allows the model to represent the changing nonlinear response of the amplifier at different excitation levels.

Values outside the measured range can either be limited to the available range or handled using extrapolation, depending on the selected configuration.

## Processing Modes

### Static processing

The function `ph_process_interp.m` estimates the level from the peak value of the complete input signal, unless a fixed level is specified manually.

One interpolated kernel set is then used for the entire signal.

```text
Complete input signal
        │
        ▼
Estimate one signal level
        │
        ▼
Select or interpolate one kernel set
        │
        ▼
Process complete signal
```

This mode is suitable for test signals and signals with a relatively stable amplitude.

### Time-varying processing

The function `time_varying_ph_process.m` divides the input into overlapping frames.

For every frame, the program:

1. estimates the current input level,
2. selects or interpolates the corresponding kernels,
3. processes all nonlinear branches,
4. applies level-dependent gain calibration,
5. reconstructs the output using overlap-add.

```text
Input signal
     │
     ▼
Overlapping frames
     │
     ▼
Level estimation for each frame
     │
     ▼
Level-dependent kernel interpolation
     │
     ▼
Parallel Hammerstein processing
     │
     ▼
Windowed overlap-add
     │
     ▼
Output signal
```

The default frame length is 64 ms, with a default hop size equal to half of the frame length.

## Oversampling and Aliasing Reduction

Polynomial nonlinearities generate new harmonics. Some of these harmonics can exceed the Nyquist frequency and produce aliasing.

To reduce this effect, nonlinear processing is performed at an increased sample rate.

The main processing sequence is:

```text
Input
  │
  ▼
FFT zero-padding oversampling
  │
  ▼
Polynomial nonlinearity
  │
  ▼
Resampling with anti-aliasing
  │
  ▼
Convolution with Hammerstein kernel
```

The oversampling factor can be configured using the `OSF` parameter.

## Gain Calibration

The function `calibrate_ph_gain_vs_level_multi.m` compares the RMS level of the model output with the RMS level of the corresponding measured amplifier response.

For each excitation level, a gain coefficient is calculated:

```text
gain = RMS measured output / RMS model output
```

These coefficients are stored in the model as a level-dependent gain curve.

During processing, the appropriate gain value is interpolated according to the currently selected input level.

This compensates for level differences between the raw Parallel Hammerstein model and the measured amplifier output.

## Repository Contents

```text
digital-tube-amplifier-model/
├── Script_PH.m
├── build.m
├── build_ph_level_model.m
├── calibrate_ph_gain_vs_level_multi.m
├── os_zeropad_fft.m
├── ph_get_kernels_for_level.m
├── ph_process_interp.m
├── plot_amp_model.m
├── time_varying_ph_process.m
├── time_varying_ph_process2.m
├── Popis Skriptu.txt
├── H3_-19dBFS.mat
├── H3_-20dBFS.mat
├── Sweep -19dBFS.wav
├── Sweep -20dBFS.wav
├── Response M3 -19dBFS.wav
└── Response M3 -20dBFS.wav
```

## File Description

### `build_ph_level_model.m`

Builds a level-dependent Parallel Hammerstein model from kernel matrices measured at multiple input levels.

The resulting model contains:

* measurement levels in dBFS,
* Hammerstein kernels,
* nonlinear branch orders,
* sampling frequency,
* optional gain-calibration data.

### `ph_get_kernels_for_level.m`

Returns the kernel matrix corresponding to a specified input level.

If the exact level is unavailable, the function interpolates all kernels between neighbouring measurement points.

### `ph_process_interp.m`

Processes a complete signal using one kernel set selected according to a fixed or estimated input level.

The function supports:

* explicit or automatically estimated input level,
* configurable oversampling,
* level-dependent gain calibration,
* optional output peak normalization.

### `time_varying_ph_process.m`

Processes the input signal frame by frame.

For each frame, it estimates the level and updates the Parallel Hammerstein kernels. Hann windows and overlap-add reconstruction are used to combine the processed frames.

### `calibrate_ph_gain_vs_level_multi.m`

Calculates gain-correction coefficients by comparing model outputs with measured amplifier responses at several excitation levels.

### `os_zeropad_fft.m`

Performs oversampling using frequency-domain zero-padding.

### `plot_amp_model.m`

Compares the measured amplifier output with the model output and provides visual analysis of their behaviour.

### `Script_PH.m`

Example processing script.

It:

* loads the generated amplifier model,
* reads an input WAV file,
* applies static or time-varying processing,
* reports input and output levels,
* plays the result,
* saves the model output to a WAV file.

## Requirements

* MATLAB
* Signal Processing Toolbox
* measured Hammerstein kernel files
* input and measured response WAV files
* a generated level-dependent model structure

## Basic Usage

### 1. Clone the repository

```bash
git clone https://github.com/SeeYaasoon/digital-tube-amplifier-model.git
cd digital-tube-amplifier-model
```

### 2. Open the project folder in MATLAB

Set the cloned repository as the current MATLAB working directory.

### 3. Prepare the amplifier model

The example script expects a generated model file containing a variable named `model`.

The model can be built using:

```matlab
model = build_ph_level_model( ...
    levelsDb, ...
    kernelFiles, ...
    referenceFrequency, ...
    sampleRate);
```

The exact kernel filenames and measurement levels must correspond to the available measurement data.

### 4. Optionally calibrate gain

```matlab
model = calibrate_ph_gain_vs_level_multi(model, levelsDb, 4);
```

### 5. Save the model

```matlab
save('amplifier_model.mat', 'model');
```

### 6. Process a signal

Update the model filename and input WAV filename in `Script_PH.m`, then run:

```matlab
Script_PH
```

The processed signal is played and written to:

```text
audio.wav
```

## Static Processing Example

```matlab
[y, levelUsed, gainApplied, outputPeak] = ph_process_interp( ...
    x, ...
    model, ...
    'OSF', 4, ...
    'NormPeak', []);
```

A fixed level can also be specified:

```matlab
[y, levelUsed] = ph_process_interp( ...
    x, ...
    model, ...
    'LevelDb', -20, ...
    'OSF', 4);
```

## Time-Varying Processing Example

```matlab
[y, levelUsed, gainApplied, outputPeak] = time_varying_ph_process( ...
    x, ...
    model, ...
    'FrameLenMs', 64, ...
    'HopMs', 32, ...
    'OSF', 4, ...
    'NormPeak', []);
```

## Measurement Data

The repository includes example excitation and response signals at selected levels:

* exponential sine sweep input signals,
* measured amplifier response signals,
* extracted Hammerstein kernel data.

These files demonstrate the expected measurement-data format, but they do not represent a complete measurement set for every possible level.

## Technologies and Topics

* MATLAB
* digital signal processing
* nonlinear system identification
* guitar amplifier modelling
* Parallel Hammerstein models
* exponential sine sweep measurements
* convolution
* polynomial nonlinearities
* impulse-response processing
* interpolation
* oversampling
* anti-aliasing
* short-time signal processing
* overlap-add reconstruction
* RMS and peak-level analysis

## Current Limitations

* the example script depends on a separately generated model `.mat` file;
* only a limited subset of measurement data is currently included in the repository;
* several filenames are configured directly inside the scripts;
* the current time-varying implementation estimates frame level using the peak value;
* abrupt changes between estimated frame levels may produce model variation;
* per-frame peak normalization, when enabled, can alter the natural dynamics of the amplifier;
* the code is intended for offline MATLAB processing rather than real-time plugin operation;
* automated tests are not currently included.

## Possible Improvements

* use RMS or a smoothed envelope instead of instantaneous frame peak level;
* add attack and release smoothing to level estimation;
* introduce hysteresis between neighbouring model levels;
* smooth interpolated kernels or model transitions between frames;
* replace per-frame normalization with physically motivated output calibration;
* add automated evaluation metrics;
* provide a complete reproducible measurement dataset;
* reorganize scripts, measurement data and generated models into separate folders;
* add plots comparing harmonic spectra and transfer characteristics;
* convert the model into a real-time C++ audio plugin.

## Project Purpose

This project was created to investigate practical methods for modelling a nonlinear tube guitar amplifier using measured input-output data.

The main objectives were:

* identifying nonlinear amplifier behaviour,
* extracting and using Hammerstein kernels,
* creating a level-dependent amplifier model,
* reducing nonlinear aliasing through oversampling,
* comparing measured and simulated responses,
* evaluating static and time-varying modelling approaches.

## Academic Context

Bachelor’s thesis project focused on the digital modelling of a nonlinear tube guitar amplifier in MATLAB.

## Author

SeeYaasoon
