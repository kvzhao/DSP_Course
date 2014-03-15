
/* Start reading here */
#include <fftw3.h>
#include "gnuplot_i.h"
#define NUM_POINTS 1024

/* Never mind this bit */
#include <stdio.h>
#include <math.h>

#define REAL 0
#define IMAG 1
#define SLEEP_LGTH  8
#define TOL     0.000001

double h_linspan[NUM_POINTS];
double freq_linspan[NUM_POINTS];

void Hf(unsigned int N, double f[NUM_POINTS], double H[NUM_POINTS]) {
    int i;
    for (i =0; i < NUM_POINTS; i++) {
        if (sin(M_PI*f[i]) < TOL )
            continue;
        H[i] = sin(M_PI*(N+1)*f[i]) / sin(M_PI*f[i]);
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

/* Resume reading here */

int main(int argc, char *argv[]) {

    int N;
    if (argc >1) {
        N = atoi(argv[1]);
    } else {
        N = 3;
    }
    gnuplot_ctrl *h;

    /* Initialize the gnuplot handle */
    printf("*** DSP example of FFT gnuplot through C ***\n") ;
    h = gnuplot_init() ;


    /* Spectrum */
    double freq_range[2] = {-3.5, 3.5};

    freq_init(freq_range, freq_linspan);
    Hf(N, freq_linspan, h_linspan);

    /* Plot the signal */
    gnuplot_resetplot(h) ;
    gnuplot_setstyle(h, "lines") ;
    printf("\n\n*** Original Signal ***\n") ;
    gnuplot_plot_xy(h,
                    freq_linspan,
                    h_linspan,
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
