# Wavelet Ping

Experimental wavelet-based ping and echo detection.

This project implements a C-based system for generating wavelets, playing them through a sound card, recording the resulting audio, and analyzing the buffer for echoes using cross-correlation and iterative signal subtraction.

## Project Structure

- `wavelet_ping_1.c`: Basic C implementation of wavelet ping and echo detection.
- `wavelet_ping_avx_2.c`: AVX-optimized implementation for faster signal subtraction.
- `mock_alsa.c`: Mock hardware layer simulating the ALSA sound library.
- `alsa/asoundlib.h`: Header for the mock ALSA library.
- `plot_results.py`: Python script for generating visual plots of the signal processing results.
- `*.png`: Visual verification graphs showing the recorded signal, initial cross-correlation peaks, and the final signal after echo removal.

## Features

- **Wavelet Generation**: Creates sine-wave wavelets of specified frequency and length.
- **Echo Detection**: Uses cross-correlation to find precise time positions and amplitudes of signal reflections.
- **Signal Subtraction**: Iteratively removes detected echoes from the recording to find subsequent, potentially overlapping reflections.
- **AVX Optimization**: Accelerates the subtraction process using SIMD instructions.
- **Hardware Simulation**: The mock ALSA layer simulates a room environment with:
  - **Direct Path**: 10ms delay.
  - **Echo 1**: 100ms delay with 50% attenuation.
  - **Echo 2**: 200ms delay with 25% attenuation.
  - **Additive Noise**: 1% random noise for robust testing.

## Compilation and Execution

### Standard Version
```bash
gcc wavelet_ping_1.c mock_alsa.c -I. -o wavelet_ping_1 -lm
./wavelet_ping_1
```

### AVX-Optimized Version
```bash
gcc wavelet_ping_avx_2.c mock_alsa.c -I. -mavx2 -o wavelet_ping_avx_2 -lm
./wavelet_ping_avx_2
```

## Visualization

To generate the test graphs, ensure you have `numpy` and `matplotlib` installed, then run the AVX version followed by the plotting script:

```bash
python3 plot_results.py
```

- `recorded_buffer.png`: Shows the raw recording with simulated echoes and noise.
- `initial_xcorr.png`: Shows the cross-correlation peaks used for echo detection.
- `final_buffer.png`: Shows the recording after all detected signals have been subtracted, effectively leaving only the noise floor.

## Dependencies

- C compiler (GCC recommended)
- Python 3 (for plotting)
- NumPy and Matplotlib (for plotting)
