
/* Start reading here */
#include <fftw3.h>
#include "gnuplot_i.h"
#define NUM_POINTS 8192

/* Never mind this bit */
#include <stdio.h>
#include <math.h>

#define SLEEP_LGTH  3
// Tolerance
#define TOL     0.0001
#define A           1

double f0, T;

double HF[NUM_POINTS];
double ht[NUM_POINTS];
double freq_linspan[NUM_POINTS];
double time_linspan[NUM_POINTS];

double Q(double f) {
    double x = 2*M_PI*T*f;
    return sin(x)/(x);
}


void gen_ht(double T, double f0, double t[NUM_POINTS], double h[NUM_POINTS]) {
    int i;
    for (i=0; i < NUM_POINTS; ++i) {
        if (t[i] > T || t[i] < -T) {
            h[i] = 0;
        } else {
            h[i] = A*cos(2*M_PI*f0*t[i]);
        }
    }
}

void gen_Hf(double f0, double f[NUM_POINTS], double H[NUM_POINTS]) {
    int i;
    for (i =0; i < NUM_POINTS; i++) {
        H[i] = A*A*T*(Q(f[i]+f0) + Q(f[i]-f0));
    }
}

void freq_init(double freq_range[2], double out[NUM_POINTS]) {
    int i;
    double f = freq_range[0];
    double f_increment = (freq_range[1]-freq_range[0]) / NUM_POINTS;
    for (i =0; i < NUM_POINTS; i++, f += f_increment) {
        out[i] = f;
    }
}

void time_init(double time_range[2], double out[NUM_POINTS]) {
    int i;
    double t = time_range[0];
    double t_increment = (time_range[1]-time_range[0]) / NUM_POINTS;
    for (i =0; i < NUM_POINTS; i++, t += t_increment) {
        out[i] = t;
    }
}
/* Resume reading here */

int main(int argc, char *argv[]) {

    if (argc >1) {
        f0 = atoi(argv[1]);
        T  = atoi(argv[2]);
    } else {
        f0 = 2;
        T  = 5;
    }
    gnuplot_ctrl *h;

    /* Initialize the gnuplot handle */
    printf("*** Digital Signal Processing through C ***\n") ;
    h = gnuplot_init() ;

    /* Spectrum */
    double freq_range[2] = {-3*f0, 3*f0};
    double time_range[2] = {-1.5*T, 1.5*T};

    freq_init(freq_range, freq_linspan);
    time_init(time_range, time_linspan);

    gen_ht(T, f0, time_linspan, ht);
    gen_Hf(f0, freq_linspan, HF);

    /* Plot the signal */
    gnuplot_resetplot(h) ;
    gnuplot_setstyle(h, "lines") ;
    gnuplot_set_xlabel(h, "Time");
    gnuplot_set_ylabel(h, "Amplitude");
    gnuplot_plot_xy(h,
                    time_linspan,
                    ht,
                    NUM_POINTS,
                    "Time Series") ;
    sleep(SLEEP_LGTH) ;

    gnuplot_resetplot(h) ;

    /* Plot the signal */
    gnuplot_resetplot(h) ;
    gnuplot_setstyle(h, "lines") ;
    gnuplot_set_xlabel(h, "Frequency");
    gnuplot_set_ylabel(h, "Amplitude");
    gnuplot_plot_xy(h,
                    freq_linspan,
                    HF,
                    NUM_POINTS,
                    "Spectrum") ;
    sleep(SLEEP_LGTH) ;

    printf("\n\n") ;
    printf("*** end of DSP example ***\n") ;
    /* Kill plotting Handler */
    gnuplot_resetplot(h);
    gnuplot_close(h);
    return 0;
}
