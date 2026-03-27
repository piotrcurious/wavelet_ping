#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>
#include "alsa/asoundlib.h"

#define SAMPLE_RATE 44100
#define PI 3.14159265
#define WAVELET_FREQ 440 // Hz
#define WAVELET_LEN 0.1 // seconds
#define BUFFER_LEN 1 // seconds
#define THRESHOLD 0.1 // minimum amplitude to detect wavelet
#define AVX_LEN 8

// Generate a wavelet of given frequency and length
float* generate_wavelet(double freq, double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    float* wavelet = malloc(n * sizeof(float)); // allocate memory
    if (!wavelet) return NULL;
    for (int i = 0; i < n; i++) {
        wavelet[i] = (float)sin(2 * PI * freq * i / SAMPLE_RATE);
    }
    return wavelet;
}

// Play a wavelet using ALSA sound library
void play_wavelet(float* wavelet, double len) {
    int n = (int) (len * SAMPLE_RATE);
    snd_pcm_t *handle;
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_PLAYBACK, 0) < 0) return;
    snd_pcm_set_params(handle, SND_PCM_FORMAT_FLOAT, SND_PCM_ACCESS_RW_INTERLEAVED, 1, SAMPLE_RATE, 1, 500000);
    snd_pcm_writei(handle, wavelet, n);
    snd_pcm_close(handle);
}

// Record a buffer of given length from sound card input using ALSA
float* record_buffer(double len) {
    int n = (int) (len * SAMPLE_RATE);
    float* buffer = malloc(n * sizeof(float));
    if (!buffer) return NULL;
    snd_pcm_t *handle;
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0) < 0) {
        free(buffer);
        return NULL;
    }
    snd_pcm_set_params(handle, SND_PCM_FORMAT_FLOAT, SND_PCM_ACCESS_RW_INTERLEAVED, 1, SAMPLE_RATE, 1, 500000);
    snd_pcm_readi(handle, buffer, n);
    snd_pcm_close(handle);
    return buffer;
}

// Find the cross-correlation between two signals
float* cross_correlate(float* signal1, int len1, float* signal2, int len2, int *out_len) {
    int len = len1 + len2 - 1;
    *out_len = len;
    float* xcorr = calloc(len, sizeof(float));
    if (!xcorr) return NULL;
    for (int i = 0; i <= len1 - len2; i++) {
        for (int j = 0; j < len2; j++) {
            xcorr[i + len2 - 1] += signal1[i + j] * signal2[j];
        }
    }
    return xcorr;
}

// Subtract a scaled and shifted signal from another signal using AVX
void subtract_signal_avx(float* x, float* y, int nx, int ny, float scale, int shift) {
    __m256 vscale = _mm256_set1_ps(scale);
    int i = 0;
    for (; i <= ny - AVX_LEN; i += AVX_LEN) {
        int buffer_idx = i + shift;
        if (buffer_idx >= 0 && buffer_idx + AVX_LEN <= nx) {
             __m256 vy = _mm256_loadu_ps(y + i);
             __m256 vx = _mm256_loadu_ps(x + buffer_idx);
             vx = _mm256_sub_ps(vx, _mm256_mul_ps(vy, vscale));
             _mm256_storeu_ps(x + buffer_idx, vx);
        } else {
            for (int j = 0; j < AVX_LEN; j++) {
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

// Analyze a buffer and find all the echoes of a wavelet
void analyze_buffer(float* buffer, float* wavelet, double len) {
    int nw = (int) (len * SAMPLE_RATE);
    int nb = (int) (BUFFER_LEN * SAMPLE_RATE);
    int count = 0;

    printf("Analyzing buffer...\n");
    while (count < 20) {
        int nr;
        float* xcorr = cross_correlate(buffer, nb, wavelet, nw, &nr);
        if (!xcorr) break;

        int max_index = -1;
        float max_value = 0;
        for (int i = 0; i < nr; i++) {
            if (xcorr[i] > max_value) {
                max_value = xcorr[i];
                max_index = i;
            }
        }

        if (max_index != -1 && max_value > THRESHOLD * (nw / 2.0)) {
            float amp = max_value / (nw / 2.0);
            int pos = max_index - nw + 1;
            printf("Echo %d: time = %f s, amplitude = %f\n", count + 1, pos / (double) SAMPLE_RATE, amp);
            subtract_signal_avx(buffer, wavelet, nb, nw, amp, pos);
            free(xcorr);
            count++;
        } else {
            free(xcorr);
            break;
        }
    }
    printf("No more echoes found. Total: %d\n", count);
}

int main() {
    float* wavelet = generate_wavelet(WAVELET_FREQ, WAVELET_LEN);
    if (!wavelet) return 1;
    play_wavelet(wavelet, WAVELET_LEN);
    float* buffer = record_buffer(BUFFER_LEN);
    if (buffer) {
        analyze_buffer(buffer, wavelet, WAVELET_LEN);
        free(buffer);
    }
    free(wavelet);
    return 0;
}
