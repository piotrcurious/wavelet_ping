#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <alsa/asoundlib.h>

#define SAMPLE_RATE 44100 // sampling rate in Hz
#define WAVELET_FREQ 440 // wavelet frequency in Hz
#define WAVELET_LEN 0.1 // wavelet length in seconds
#define BUFFER_LEN 1 // buffer length in seconds
#define THRESHOLD 0.01 // threshold for echo detection
#define AVX_LEN 8 // number of floats processed by AVX

// a function to generate a wavelet of a given frequency and length
float* generate_wavelet(float freq, float len) {
    int n = (int) (len * SAMPLE_RATE); // number of samples
    float* wavelet = (float*) malloc(n * sizeof(float)); // allocate memory
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
    snd_pcm_hw_params_t* params; // hardware parameters
    snd_pcm_uframes_t frames; // number of frames per period

    // open PCM device for playback
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_PLAYBACK, 0) < 0) {
        fprintf(stderr, "Error: unable to open PCM device\n");
        exit(1);
    }

    // allocate hardware parameters object and fill it with default values
    snd_pcm_hw_params_alloca(&params);
    snd_pcm_hw_params_any(handle, params);

    // set parameters: interleaved mode, 32-bit little-endian format, 1 channel, sampling rate
    snd_pcm_hw_params_set_access(handle, params, SND_PCM_ACCESS_RW_INTERLEAVED);
    snd_pcm_hw_params_set_format(handle, params, SND_PCM_FORMAT_FLOAT_LE);
    snd_pcm_hw_params_set_channels(handle, params, 1);
    snd_pcm_hw_params_set_rate(handle, params, SAMPLE_RATE, 0);

    // write parameters to the device and prepare it for use
    if (snd_pcm_hw_params(handle, params) < 0) {
        fprintf(stderr, "Error: unable to set hardware parameters\n");
        exit(1);
    }

    // get the number of frames per period and the period time
    snd_pcm_hw_params_get_period_size(params, &frames, NULL);
    snd_pcm_hw_params_get_period_time(params, NULL, NULL);

    int n = (int) (len * SAMPLE_RATE); // number of samples to play
    int i = 0; // current sample index

    while (i < n) {
        // write frames to the device buffer
        int err = snd_pcm_writei(handle, wavelet + i, frames);
        if (err == -EPIPE) {
            // buffer underrun occurred, try to recover
            snd_pcm_prepare(handle);
        } else if (err < 0) {
            fprintf(stderr, "Error: write failed\n");
            exit(1);
        } else {
            i += err; // update sample index
        }
    }

    // wait for the device to finish playing and close it
    snd_pcm_drain(handle);
    snd_pcm_close(handle);
}

// a function to record a buffer using ALSA
float* record_buffer(float len) {
    snd_pcm_t* handle; // PCM device handle
    snd_pcm_hw_params_t* params; // hardware parameters
    snd_pcm_uframes_t frames; // number of frames per period

    // open PCM device for capture
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0) < 0) {
        fprintf(stderr, "Error: unable to open PCM device\n");
        exit(1);
    }

    // allocate hardware parameters object and fill it with default values
    snd_pcm_hw_params_alloca(&params);
    snd_pcm_hw_params_any(handle, params);

    // set parameters: interleaved mode, 32-bit little-endian format, 1 channel, sampling rate
    snd_pcm_hw_params_set_access(handle, params, SND_PCM_ACCESS_RW_INTERLEAVED);
    snd_pcm_hw_params_set_format(handle, params, SND_PCM_FORMAT_FLOAT_LE);
    snd_pcm_hw_params_set_channels(handle, params, 1);
    snd_pcm_hw_params_set_rate(handle, params, SAMPLE_RATE, 0);

    // write parameters to the device and prepare it for use
    if (snd_pcm_hw_params(handle, params) < 0) {
        fprintf(stderr, "Error: unable to set hardware parameters\n");
        exit(1);
    }

    // get the number of frames per period and the period time
    snd_pcm_hw_params_get_period_size(params, &frames, NULL);
    snd_pcm_hw_params_get_period_time(params, NULL, NULL);

    int n = (int) (len * SAMPLE_RATE); // number of samples to record
    float* buffer = (float*) malloc(n * sizeof(float)); // allocate memory
    int i = 0; // current sample index

    while (i < n) {
        // read frames from the device buffer
        int err = snd_pcm_readi(handle, buffer + i, frames);
        if (err == -EPIPE) {
            // buffer overrun occurred, try to recover
            snd_pcm_prepare(handle);
        } else if (err < 0) {
            fprintf(stderr, "Error: read failed\n");
            exit(1);
        } else {
            i += err; // update sample index
        }
    }

    // close the device
    snd_pcm_close(handle);

    return buffer;
}

// a function to compute the cross-correlation between two signals
float* cross_correlation(float* x, float* y, int nx, int ny) {
    int n = nx + ny - 1; // length of the cross-correlation
    float* r = (float*) malloc(n * sizeof(float)); // allocate memory
    for (int k = 0; k < n; k++) {
        r[k] = 0; // initialize to zero
        int i_min = (k >= ny - 1) ? k - (ny - 1) : 0; // lower bound of i
        int i_max = (k < nx - 1) ? k : nx - 1; // upper bound of i
        for (int i = i_min; i <= i_max; i++) {
            r[k] += x[i] * y[k - i]; // sum the products
        }
    }
    return r;
}

// a function to find the maximum value and its index in an array
void find_max(float* x, int n, float* max_val, int* max_idx) {
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
    __m256 vx, vy, vscale; // AVX registers for x, y and scale
    vscale = _mm256_set1_ps(scale); // broadcast scale to all elements of vscale
    for (int i = 0; i < nx; i += AVX_LEN) {
        vx = _mm256_loadu_ps(x + i); // load AVX_LEN elements of x from memory
        if (i + shift >= 0 && i + shift + AVX_LEN <= ny) {
            vy = _mm256_loadu_ps(y + i + shift); // load AVX_LEN elements of y from memory
        } else {
            float temp[AVX_LEN]; // temporary array to store y elements
            for (int j = 0; j < AVX_LEN; j++) {
                if (i + shift + j >= 0 && i + shift + j < ny) {
                    temp[j] = y[i + shift + j]; // copy valid y element
                } else {
                    temp[j] = 0; // pad with zero
                }
            }
            vy = _mm256_loadu_ps(temp); // load temp array to AVX register
        }
        vy = _mm256_mul_ps(vy, vscale); // multiply y by scale
        vx = _mm256_sub_ps(vx, vy); // subtract y from x
        _mm256_storeu_ps(x + i, vx); // store the result back to x in memory
    }
}

// a function to analyze a buffer and find the echoes of a wavelet using cross-correlation and subtraction
void analyze_buffer(float* buffer, float* wavelet, float buffer_len, float wavelet_len) {
    int nb = (int) (buffer_len * SAMPLE_RATE); // number of samples in buffer
    int nw = (int) (wavelet_len * SAMPLE_RATE); // number of samples in wavelet

    printf("Analyzing buffer...\n

          // compute the cross-correlation between buffer and wavelet
    float* r = cross_correlation(buffer, wavelet, nb, nw);

    // find the maximum value and its index in the cross-correlation
    float max_val;
    int max_idx;
    find_max(r, nb + nw - 1, &max_val, &max_idx);

    // loop until no more echoes are found
    int count = 0; // number of echoes found
    while (max_val > THRESHOLD) {
        // calculate the time and amplitude of the echo
        float time = (float) (max_idx - nw + 1) / SAMPLE_RATE;
        float amp = max_val / (nw * nw);

        // print the echo information
        printf("Echo %d: time = %.3f s, amplitude = %.3f\n", count + 1, time, amp);

        // subtract the scaled and shifted wavelet from the buffer
        subtract_signal(buffer, wavelet, nb, nw, amp, max_idx - nw + 1);

        // compute the cross-correlation again
        free(r); // free the previous cross-correlation memory
        r = cross_correlation(buffer, wavelet, nb, nw);

        // find the maximum value and its index in the cross-correlation
        find_max(r, nb + nw - 1, &max_val, &max_idx);

        // increment the echo count
        count++;
    }

    // free the memory allocated for cross-correlation
    free(r);

    printf("No more echoes found. Total number of echoes: %d\n", count);
}

// the main function
int main() {
    // generate a wavelet of a given frequency and length
    float* wavelet = generate_wavelet(WAVELET_FREQ, WAVELET_LEN);

    // play the wavelet using ALSA
    play_wavelet(wavelet, WAVELET_LEN);

    // record a buffer using ALSA
    float* buffer = record_buffer(BUFFER_LEN);

    // analyze the buffer and find the echoes of the wavelet
    analyze_buffer(buffer, wavelet, BUFFER_LEN, WAVELET_LEN);

    // free the memory allocated for wavelet and buffer
    free(wavelet);
    free(buffer);

    return 0;
}
