#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "alsa/asoundlib.h"

#define SAMPLE_RATE 44100
#define PI 3.14159265
#define WAVELET_FREQ 1000 // Hz
#define WAVELET_LEN 0.1 // seconds
#define BUFFER_LEN 1 // seconds
#define THRESHOLD 0.1 // minimum amplitude to detect wavelet

// Generate a wavelet of given frequency and length
float* generate_wavelet(double freq, double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    float* wavelet = malloc(n * sizeof(float)); // allocate memory
    for (int i = 0; i < n; i++) {
        // use sine function to generate wavelet
        wavelet[i] = (float)sin(2 * PI * freq * i / SAMPLE_RATE);
    }
    return wavelet;
}

// Play a wavelet using ALSA sound library
void play_wavelet(float* wavelet, double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    int err;
    snd_pcm_t *handle;

    // open PCM device for playback
    err = snd_pcm_open(&handle, "default", SND_PCM_STREAM_PLAYBACK, 0);
    if (err < 0) {
        fprintf(stderr, "Unable to open PCM device: %s\n", snd_strerror(err));
        exit(1);
    }

    // set hardware parameters
    err = snd_pcm_set_params(handle,
                             SND_PCM_FORMAT_FLOAT,
                             SND_PCM_ACCESS_RW_INTERLEAVED,
                             1, // number of channels
                             SAMPLE_RATE,
                             1, // allow software resampling
                             500000); // latency in microseconds
    if (err < 0) {
        fprintf(stderr, "Unable to set PCM parameters: %s\n", snd_strerror(err));
        exit(1);
    }

    // write wavelet data to PCM device and play it
    snd_pcm_writei(handle, wavelet, n);

    // close PCM device
    snd_pcm_close(handle);
}

// Record a buffer of given length from sound card input using ALSA sound library
float* record_buffer(double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    float* buffer = malloc(n * sizeof(float)); // allocate memory
    int err;
    snd_pcm_t *handle;

    // open PCM device for capture
    err = snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0);
    if (err < 0) {
        fprintf(stderr, "Unable to open PCM device: %s\n", snd_strerror(err));
        exit(1);
    }

    // set hardware parameters
    err = snd_pcm_set_params(handle,
                             SND_PCM_FORMAT_FLOAT,
                             SND_PCM_ACCESS_RW_INTERLEAVED,
                             1, // number of channels
                             SAMPLE_RATE,
                             1, // allow software resampling
                             500000); // latency in microseconds
    if (err < 0) {
        fprintf(stderr, "Unable to set PCM parameters: %s\n", snd_strerror(err));
        exit(1);
    }

    // read buffer data from PCM device and record it
    snd_pcm_readi(handle, buffer, n);

    // close PCM device
    snd_pcm_close(handle);

    return buffer;
}

// Find the cross-correlation between two signals
float* cross_correlate(float* signal1, int len1, float* signal2, int len2, int *out_len) {
    // allocate memory for cross-correlation array
    int len = len1 + len2 - 1; // length of cross-correlation array
    *out_len = len;
    float* xcorr = calloc(len, sizeof(float));

    for (int i = 0; i <= len1 - len2; i++) {
        for (int j = 0; j < len2; j++) {
            xcorr[i + len2 - 1] += signal1[i + j] * signal2[j];
        }
    }
    return xcorr;
}

// Find the time position of a wavelet in a buffer using cross-correlation and thresholding
int find_wavelet(float* buffer, float* wavelet, double len, float *out_amp) {
    int n = (int) (len * SAMPLE_RATE); // number of samples in wavelet
    int m = (int) (BUFFER_LEN * SAMPLE_RATE); // number of samples in buffer
    int xcorr_len;
    float* xcorr = cross_correlate(buffer, m, wavelet, n, &xcorr_len); // cross-correlate buffer and wavelet
    int max_index = -1; // index of maximum value in xcorr
    float max_value = 0; // maximum value in xcorr
    for (int i = 0; i < xcorr_len; i++) {
        // find the maximum value and its index in xcorr
        if (xcorr[i] > max_value) {
            max_value = xcorr[i];
            max_index = i;
        }
    }

    if (max_value > THRESHOLD * (n / 2.0)) { // Normalize threshold by wavelet energy
        *out_amp = max_value / (n / 2.0);
        free(xcorr); // free memory allocated for xcorr
        return max_index - n + 1;
    } else {
        free(xcorr); // free memory allocated for xcorr
        return -1;
    }
}

// Subtract a wavelet from a buffer at a given time position
void subtract_wavelet(float* buffer, float* wavelet, double len, int pos, float amp) {
    int n = (int) (len * SAMPLE_RATE); // number of samples in wavelet
    int m = (int) (BUFFER_LEN * SAMPLE_RATE);
    for (int i = 0; i < n; i++) {
        if (pos + i >= 0 && pos + i < m) {
            // subtract the wavelet from the buffer element-wise at the given position
            buffer[pos + i] -= amp * wavelet[i];
        }
    }
}

// Analyze a buffer and find all the echoes of a wavelet using subtraction and repetition
void analyze_buffer(float* buffer, float* wavelet, double len) {
    float amp;
    int pos = find_wavelet(buffer, wavelet, len, &amp); // find the first wavelet in the buffer
    int count = 0;
    while (pos != -1 && count < 20) {
        // while there is a wavelet in the buffer, do the following:
        printf("Found wavelet at %f seconds, amp: %f\n", pos / (double) SAMPLE_RATE, amp); // print the time position of the wavelet
        subtract_wavelet(buffer, wavelet, len, pos, amp); // subtract the wavelet from the buffer
        pos = find_wavelet(buffer, wavelet, len, &amp); // find the next wavelet in the buffer
        count++;
    }
}

// Main function
int main() {
    float* wavelet = generate_wavelet(WAVELET_FREQ, WAVELET_LEN); // generate a wavelet of given frequency and length
    if (!wavelet) return 1;
    play_wavelet(wavelet, WAVELET_LEN); // play the wavelet using sound card output
    float* buffer = record_buffer(BUFFER_LEN); // record a buffer of given length from sound card input
    if (buffer) {
        analyze_buffer(buffer, wavelet, WAVELET_LEN); // analyze the buffer and find all the echoes of the wavelet
        free(buffer); // free memory allocated for buffer
    }

    free(wavelet); // free memory allocated for wavelet

    return 0;
}
