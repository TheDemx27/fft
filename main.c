#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>

void die(const char *message)
{
    if(errno) {
        perror(message);
    } else {
        printf("ERROR: %s\n", message);
    }
    exit(1);
}

/*Window Functions*/

typedef void (*cb_func)(int n, int length, float *window);

void genwindow(int length, float *window, cb_func winfunc)
{
    for(int n = 0; n < length; n++){
        winfunc(n, length, window);
    }
}

void boxcar(int n, int length, float *window)
{
    window[n] = 1;
}

void hamming(int n, int length, float *window)
{
    window[n] = 0.54 - 0.46*cos(2*M_PI*n/length);
}

void hanning(int n, int length, float *window)
{
    window[n] = 0.5 - 0.5*cos(2*M_PI*n/length);
}

void blackman(int n, int length, float *window)
{
    window[n] = 0.42 + 0.5*cos(2*M_PI*n/(length - 1)) + 0.08*cos(4*M_PI*n/(length - 1));
}

/*End Window Functions*/

float *gensample(int sampleLength, float samplePrec, cb_func win)
{
    float *window_samples = malloc(sizeof(float)*sampleLength);
    float *samples = malloc(sizeof(float)*sampleLength);
    if (samples == NULL) die("MALLOC ERROR");

    genwindow(sampleLength, window_samples, hamming); // edits window_samples

    for(int i = 0; i < sampleLength; i++) {
        samples[i] = window_samples[i]*(sin(30.5*2*M_PI*i*1/samplePrec)+sin(15.5*2*M_PI*i*1/samplePrec));
    }
    return samples;
}

float *dft(float sample[], float length, float *X_real, float *X_img, float *K)
{
    float *result = malloc(length*sizeof(float));
    if (result == NULL) die("MALLOC ERROR");
    float x_real;
    float x_img;
    float arg;
    float k_f = -length/2;

    for(int k = 0; k < length; k++) { // freq
        for(int n = 0; n < length; n++) { // time
            arg = -2*M_PI*k_f*n/length;

            x_real += sample[n]*cos(arg);
            x_img += sample[n]*sin(arg);
        }

        X_real[k] = x_real;
        X_img[k] = x_img;
        K[k] = k_f;

        result[k] = sqrt(pow(x_real,2)+pow(x_img,2));

        k_f += 1;
        x_real = 0;
        x_img = 0;
    }
    return result;
}

void log_data(int length, float *result, float *K)
{
    FILE *fp;
    fp = fopen("data.txt", "w");
    for(int i = 0; i < length; i++) {
        fprintf(fp, "%f, %f, ", K[i], result[i]);
    }
    fclose(fp);
}

int main(int argc, char *argv[])
{
    if(argc < 3) die("Insufficient arguments supplied.");
    int length = atoll(argv[1]);
    float prec = atof(argv[2]);

    printf("sizeof(length) = %lu\n", sizeof(length));

    float *samples = gensample(length, prec, hamming);
    float *K = malloc(length*sizeof(float));
    memset(K, 0, length);
    float *X_real = malloc(length*sizeof(float));
    memset(X_real, 0, length);
    float *X_img = malloc(length*sizeof(float));
    memset(X_img, 0, length);

    float *result = dft(samples, length, X_real, X_img, K);
    log_data(length, result, K);

    free(X_real);
    free(X_img);
    free(K);
    free(samples);
    return 0;
}
