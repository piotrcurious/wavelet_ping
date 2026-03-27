#include "alsa/asoundlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SAMPLES (44100 * 5) // 5 seconds
float global_playback_buffer[MAX_SAMPLES];
int playback_len = 0;

struct _snd_pcm {
    int stream;
};

int snd_pcm_open(snd_pcm_t **pcm, const char *name, snd_pcm_stream_t stream, int mode) {
    *pcm = (snd_pcm_t *)malloc(sizeof(snd_pcm_t));
    (*pcm)->stream = stream;
    if (stream == SND_PCM_STREAM_PLAYBACK) {
        playback_len = 0;
        memset(global_playback_buffer, 0, sizeof(global_playback_buffer));
    }
    return 0;
}

int snd_pcm_close(snd_pcm_t *pcm) {
    free(pcm);
    return 0;
}

const char *snd_strerror(int errnum) {
    return "ALSA error";
}

int snd_pcm_set_params(snd_pcm_t *pcm, snd_pcm_format_t format, snd_pcm_access_t access, unsigned int channels, unsigned int rate, int soft_resample, unsigned int latency) {
    return 0;
}

snd_pcm_sframes_t snd_pcm_writei(snd_pcm_t *pcm, const void *buffer, snd_pcm_uframes_t size) {
    int to_copy = (size < MAX_SAMPLES - playback_len) ? size : (MAX_SAMPLES - playback_len);
    memcpy(global_playback_buffer + playback_len, buffer, to_copy * sizeof(float));
    playback_len += to_copy;
    return to_copy;
}

snd_pcm_sframes_t snd_pcm_readi(snd_pcm_t *pcm, void *buffer, snd_pcm_uframes_t size) {
    float *out = (float *)buffer;
    for (int i = 0; i < size; i++) {
        out[i] = 0;
        // Direct sound (with a small delay)
        int delay0 = 441; // 10ms
        if (i >= delay0 && i - delay0 < playback_len) {
            out[i] += global_playback_buffer[i - delay0];
        }
        // Echo 1
        int delay1 = 4410; // 100ms
        if (i >= delay1 && i - delay1 < playback_len) {
            out[i] += 0.5f * global_playback_buffer[i - delay1];
        }
        // Echo 2
        int delay2 = 8820; // 200ms
        if (i >= delay2 && i - delay2 < playback_len) {
            out[i] += 0.25f * global_playback_buffer[i - delay2];
        }

        // Add some noise
        float noise = ((float)rand() / (float)RAND_MAX - 0.5f) * 0.01f;
        out[i] += noise;
    }
    return size;
}

int snd_pcm_recover(snd_pcm_t *pcm, int err, int silent) { return 0; }
int snd_pcm_prepare(snd_pcm_t *pcm) { return 0; }
int snd_pcm_drain(snd_pcm_t *pcm) { return 0; }

int snd_pcm_hw_params_any(snd_pcm_t *pcm, snd_pcm_hw_params_t *params) { return 0; }
int snd_pcm_hw_params_set_access(snd_pcm_t *pcm, snd_pcm_hw_params_t *params, snd_pcm_access_t _access) { return 0; }
int snd_pcm_hw_params_set_format(snd_pcm_t *pcm, snd_pcm_hw_params_t *params, snd_pcm_format_t val) { return 0; }
int snd_pcm_hw_params_set_channels(snd_pcm_t *pcm, snd_pcm_hw_params_t *params, unsigned int val) { return 0; }
int snd_pcm_hw_params_set_rate(snd_pcm_t *pcm, snd_pcm_hw_params_t *params, unsigned int val, int dir) { return 0; }
int snd_pcm_hw_params(snd_pcm_t *pcm, snd_pcm_hw_params_t *params) { return 0; }
int snd_pcm_hw_params_get_period_size(const snd_pcm_hw_params_t *params, snd_pcm_uframes_t *frames, int *dir) {
    *frames = 1024;
    return 0;
}
int snd_pcm_hw_params_get_period_time(const snd_pcm_hw_params_t *params, unsigned int *val, int *dir) {
    *val = 23219;
    return 0;
}
