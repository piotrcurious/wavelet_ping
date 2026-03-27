#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>
#include "alsa/asoundlib.h"

#define SAMPLE_RATE 44100 // sampling rate in Hz
#define WAVELET_FREQ 440 // wavelet frequency in Hz
#define WAVELET_LEN 0.1 // wavelet length in seconds
#define BUFFER_LEN 1 // buffer length in seconds
#define THRESHOLD 0.1 // threshold for echo detection
#define AVX_LEN 8 // number of floats processed by AVX

// a function to generate a wavelet of a given frequency and length
float* generate_wavelet(float freq, float len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    float* wavelet = (float*) malloc(n * sizeof(float)); // allocate memory
    if (!wavelet) return NULL;
    float omega = 2 * M_PI * freq; // angular frequency
    for (int i = 0; i < n; i++) {
        float t = (float) i / SAMPLE_RATE; // time
        wavelet[i] = sin(omega * t); // sine wave
    }
    return wavelet;
}

// a function to play a wavelet using ALSA
void play_wavelet(float* wavelet, float len) {
    snd_pcm_t* handle; // PCM device handle

    // open PCM device for playback
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_PLAYBACK, 0) < 0) {
        fprintf(stderr, "Error: unable to open PCM device\n");
        exit(1);
    }

    int n = (int) (len * SAMPLE_RATE); // number of samples to play
    snd_pcm_writei(handle, wavelet, n);

    // wait for the device to finish playing and close it
    snd_pcm_drain(handle);
    snd_pcm_close(handle);
}

// a function to record a buffer using ALSA
float* record_buffer(float len) {
    snd_pcm_t* handle; // PCM device handle

    // open PCM device for capture
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0) < 0) {
        fprintf(stderr, "Error: unable to open PCM device\n");
        exit(1);
    }

    int n = (int) (len * SAMPLE_RATE); // number of samples to record
    float* buffer = (float*) malloc(n * sizeof(float)); // allocate memory
    if (!buffer) return NULL;

    snd_pcm_readi(handle, buffer, n);

    // close the device
    snd_pcm_close(handle);

    return buffer;
}

// a function to compute the cross-correlation between two signals
float* cross_correlation(float* x, float* y, int nx, int ny, int *out_len) {
    int n = nx + ny - 1; // length of the cross-correlation
    *out_len = n;
    float* r = (float*) calloc(n, sizeof(float)); // allocate memory and zero it
    if (!r) return NULL;
    for (int k = 0; k <= nx - ny; k++) {
        for (int i = 0; i < ny; i++) {
            r[k + ny - 1] += x[k + i] * y[i];
        }
    }
    return r;
}

// a function to find the maximum value and its index in an array
void find_max(float* x, int n, float* max_val, int* max_idx) {
    if (n <= 0) {
        *max_val = 0;
        *max_idx = -1;
        return;
    }
    *max_val = x[0]; // initialize max value
    *max_idx = 0; // initialize max index
    for (int i = 1; i < n; i++) {
        if (x[i] > *max_val) {
            *max_val = x[i]; // update max value
            *max_idx = i; // update max index
        }
    }
}

// a function to subtract a scaled and shifted signal from another signal using AVX
void subtract_signal(float* x, float* y, int nx, int ny, float scale, int shift) {
    __m256 vscale = _mm256_set1_ps(scale); // broadcast scale to all elements of vscale
    int i = 0;
    for (; i <= ny - AVX_LEN; i += AVX_LEN) {
        int buffer_idx = i + shift;
        if (buffer_idx >= 0 && buffer_idx + AVX_LEN <= nx) {
             __m256 vy = _mm256_loadu_ps(y + i);
             __m256 vx = _mm256_loadu_ps(x + buffer_idx);
             vx = _mm256_sub_ps(vx, _mm256_mul_ps(vy, vscale));
             _mm256_storeu_ps(x + buffer_idx, vx);
        } else {
            // Edge of buffer
            for (int j = 0; j < AVX_LEN && i + j < ny; j++) {
                if (shift + i + j >= 0 && shift + i + j < nx) {
                    x[shift + i + j] -= y[i + j] * scale;
                }
            }
        }
    }
    // Remainder
    for (; i < ny; i++) {
        if (shift + i >= 0 && shift + i < nx) {
            x[shift + i] -= y[i] * scale;
        }
    }
}

// a function to analyze a buffer and find the echoes of a wavelet using cross-correlation and subtraction
void analyze_buffer(float* buffer, float* wavelet, float buffer_len, float wavelet_len) {
    int nb = (int) (buffer_len * SAMPLE_RATE); // number of samples in buffer
    int nw = (int) (wavelet_len * SAMPLE_RATE); // number of samples in wavelet

    // Save initial buffer for plotting
    FILE* f = fopen("recorded_buffer.bin", "wb");
    if (f) {
        fwrite(buffer, sizeof(float), nb, f);
        fclose(f);
    }

    printf("Analyzing buffer...\n");

    // compute the cross-correlation between buffer and wavelet
    int nr;
    float* r = cross_correlation(buffer, wavelet, nb, nw, &nr);
    if (!r) return;

    // Save initial cross-correlation for plotting
    f = fopen("initial_xcorr.bin", "wb");
    if (f) {
        fwrite(r, sizeof(float), nr, f);
        fclose(f);
    }

    // find the maximum value and its index in the cross-correlation
    float max_val;
    int max_idx;
    find_max(r, nr, &max_val, &max_idx);

    // loop until no more echoes are found
    int count = 0; // number of echoes found
    while (max_idx != -1 && max_val > THRESHOLD * (nw / 2.0)) {
        // calculate the time and amplitude of the echo
        float time = (float) (max_idx - nw + 1) / SAMPLE_RATE;
        float amp = max_val / (nw / 2.0);

        // print the echo information
        printf("Echo %d: time = %.3f s, amplitude = %.3f\n", count + 1, time, amp);

        // subtract the scaled and shifted wavelet from the buffer
        subtract_signal(buffer, wavelet, nb, nw, amp, max_idx - nw + 1);

        // compute the cross-correlation again
        free(r); // free the previous cross-correlation memory
        r = cross_correlation(buffer, wavelet, nb, nw, &nr);
        if (!r) break;

        // find the maximum value and its index in the cross-correlation
        find_max(r, nr, &max_val, &max_idx);

        // increment the echo count
        count++;
        if (count > 20) break; // Safety break
    }

    // Save final buffer for plotting
    f = fopen("final_buffer.bin", "wb");
    if (f) {
        fwrite(buffer, sizeof(float), nb, f);
        fclose(f);
    }

    // free the memory allocated for cross-correlation
    if (r) free(r);

    printf("No more echoes found. Total number of echoes: %d\n", count);
}

// the main function
int main() {
    // generate a wavelet of a given frequency and length
    float* wavelet = generate_wavelet(WAVELET_FREQ, WAVELET_LEN);
    if (!wavelet) return 1;

    // play the wavelet using ALSA
    play_wavelet(wavelet, WAVELET_LEN);

    // record a buffer using ALSA
    float* buffer = record_buffer(BUFFER_LEN);

    // analyze the buffer and find the echoes of the wavelet
    if (buffer) {
        analyze_buffer(buffer, wavelet, BUFFER_LEN, WAVELET_LEN);
        free(buffer);
    }

    // free the memory allocated for wavelet and buffer
    free(wavelet);

    return 0;
}
