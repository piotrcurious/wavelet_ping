#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <alsa/asoundlib.h>

#define SAMPLE_RATE 44100
#define PI 3.14159265
#define WAVELET_FREQ 1000 // Hz
#define WAVELET_LEN 0.1 // seconds
#define BUFFER_LEN 1 // seconds
#define THRESHOLD 0.01 // minimum amplitude to detect wavelet

// Generate a wavelet of given frequency and length
double* generate_wavelet(double freq, double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    double* wavelet = malloc(n * sizeof(double)); // allocate memory
    for (int i = 0; i < n; i++) {
        // use sine function to generate wavelet
        wavelet[i] = sin(2 * PI * freq * i / SAMPLE_RATE);
    }
    return wavelet;
}

// Play a wavelet using ALSA sound library
void play_wavelet(double* wavelet, double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    int err;
    snd_pcm_t *handle;
    snd_pcm_sframes_t frames;

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
    frames = snd_pcm_writei(handle, wavelet, n);
    if (frames < 0) {
        frames = snd_pcm_recover(handle, frames, 0);
    }
    if (frames < 0) {
        fprintf(stderr, "Unable to write to PCM device: %s\n", snd_strerror(err));
        exit(1);
    }

    // close PCM device
    snd_pcm_close(handle);
}

// Record a buffer of given length from sound card input using ALSA sound library
double* record_buffer(double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    double* buffer = malloc(n * sizeof(double)); // allocate memory
    int err;
    snd_pcm_t *handle;
    snd_pcm_sframes_t frames;

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
    frames = snd_pcm_readi(handle, buffer, n);
    if (frames < 0) {
        frames = snd_pcm_recover(handle, frames, 0);
    }
    if (frames < 0) {
        fprintf(stderr, "Unable to read from PCM device: %s\n", snd_strerror(err));
        exit(1);
    }

    // close PCM device
    snd_pcm_close(handle);

    return buffer;
}

// Find the cross-correlation between two signals of different lengths
double* cross_correlate(double* signal1, int len1, double* signal2, int len2) {
    // allocate memory for cross-correlation array
    int len = len1 + len2 - 1; // length of cross-correlation array
    double* xcorr = malloc(len * sizeof(double));

    // compute cross-correlation by sliding signal2 over signal1 and multiplying element-wise
    for (int i = 0; i < len; i++) {
        xcorr[i] = 0; // initialize to zero
        int start = i - len2 + 1; // start index of signal1
        if (start < 0) start = 0; // clamp to zero if negative
        int end = i; // end index of signal1
        if (end >= len1) end = len1 - 1; // clamp to len1 - 1 if too large
        for (int j = start; j <= end; j++) {
            // multiply corresponding elements of signal1 and signal2 and add to xcorr[i]
            xcorr[i] += signal1[j] * signal2[i - j];
        }
    }
    return xcorr;
}

// Find the time position of a wavelet in a buffer using cross-correlation and thresholding
int find_wavelet(double* buffer, double* wavelet, double len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples in wavelet
    int m = (int) (BUFFER_LEN * SAMPLE_RATE); // number of samples in buffer
    double* xcorr = cross_correlate(buffer, m, wavelet, n); // cross-correlate buffer and wavelet
    int max_index = -1; // index of maximum value in xcorr
    double max_value = 0; // maximum value in xcorr
    for (int i = 0; i < m + n - 1; i++) {
        // find the maximum value and its index in xcorr
        if (xcorr[i] > max_value) {
            max_value = xcorr[i];
            max_index = i;
        }
    }
    free(xcorr); // free memory allocated for xcorr
    if (max_value > THRESHOLD) {
        // if the maximum value is above the threshold, return the time position of the wavelet
        return max_index - n + 1;
    } else {
        // otherwise, return -1 to indicate no wavelet found
        return -1;
    }
}

// Subtract a wavelet from a buffer at a given time position
void subtract_wavelet(double* buffer, double* wavelet, double len, int pos) {
    int n = (int) (len * SAMPLE_RATE); // number of samples in wavelet
    for (int i = 0; i < n; i++) {
        // subtract the wavelet from the buffer element-wise at the given position
        buffer[pos + i] -= wavelet[i];
    }
}

// Analyze a buffer and find all the echoes of a wavelet using subtraction and repetition
void analyze_buffer(double* buffer, double* wavelet, double len) {
    int pos = find_wavelet(buffer, wavelet, len); // find the first wavelet in the buffer
    while (pos != -1) {
        // while there is a wavelet in the buffer, do the following:
        printf("Found wavelet at %f seconds\n", pos / (double) SAMPLE_RATE); // print the time position of the wavelet
        subtract_wavelet(buffer, wavelet, len, pos); // subtract the wavelet from the buffer
        pos = find_wavelet(buffer, wavelet, len); // find the next wavelet in the buffer
    }
}

// Main function
int main() {
    double* wavelet = generate_wavelet(WAVELET_FREQ, WAVELET_LEN); // generate a wavelet of given frequency and length
    play_wavelet(wavelet, WAVELET_LEN); // play the wavelet using sound card output
    sleep(1); // wait for one second before recording
    double* buffer = record_buffer(BUFFER_LEN); // record a buffer of given length from sound card input
    analyze_buffer(buffer, wavelet, WAVELET_LEN); // analyze the buffer and find all the echoes of the wavelet

    free(wavelet); // free memory allocated for wavelet
    free(buffer); // free memory allocated for buffer

    return 0;
}

