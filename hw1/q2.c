/* Start reading here */
#include <fftw3.h>
#include "gnuplot_i.h"
#define NUM_POINTS 1024

/* Never mind this bit */
#include <stdio.h>
#include <math.h>

unsigned int EXPAN_TERM = 10;
unsigned int SLEEP_LGTH  = 3;

double time_series     [NUM_POINTS];
double square_wave     [NUM_POINTS];
double triangle_wave   [NUM_POINTS];
double sawtooth_wave   [NUM_POINTS];

void init_time_series(double t_range[2]) {
    double time_length = t_range[1] - t_range[0];
    double time_step = time_length / NUM_POINTS;
    double t;
    int i;
    for (i =0, t = t_range[0]; i< NUM_POINTS; ++i, t += time_step) {
        time_series[i] = t;
    }
}

void init_square_wave(double t_range[2],
                      double period,
                      double ts[NUM_POINTS],
                      double f[NUM_POINTS])
{
    double time_length = t_range[1] - t_range[0];
    double time_step = time_length / NUM_POINTS;
    double t;
    int i, A = 1;
    int pn = period / time_step, p2 = pn >> 2;

    for (i =0, t = t_range[0]; i< NUM_POINTS; i++, t += time_step) {
        if ((i%pn) > p2)
            f[i] = A;
        else
            f[i] = -A;
        ts[i] = t;
    }
}
/* Resume reading here */
#define SQUARE_WAVE 1
#define TRIANGLE_WAVE 2
#define SAWTOOTH_WAVE 3

void fourier_expansion(unsigned int sig_type) {
    switch (sig_type) {
        case SQUARE_WAVE:
        {
            // i: time index, k: expansion term
            int i, k;
            // linear combination of sinusoids
            double lincom =0;
            /* Whole time series */
            for (i =0; i < NUM_POINTS; ++i) {
                /* Expansion term */
                double t = time_series[i];
                for (k =0; k < EXPAN_TERM; ++k) {
                    lincom += sin(2*M_PI*(2*k+1) * t) / (2*k+1);
                }
                square_wave[i] = 4*lincom/M_PI;
                // reset the accumulation term
                lincom =0;
            }
        }
        break;

        case TRIANGLE_WAVE:
        {
            // i: time index, k: expansion term
            int i, k;
            // linear combination of sinusoids
            double lincom =0;
            /* Whole time series */
            for (i =0; i < NUM_POINTS; ++i) {
                /* Expansion term */
                double t = time_series[i];
                for (k =0; k < EXPAN_TERM; ++k) {
                    lincom += pow(-1.0,k) * sin(2*M_PI*(2*k+1) * t) /(2*k+1)/(2*k+1);
                }
                triangle_wave[i] = 4*lincom /(M_PI*M_PI);
                // reset the accumulation term
                lincom =0;
            }
        }
        break;

        case SAWTOOTH_WAVE:
        {
            // i: time index, k: expansion term
            int i, k;
            // linear combination of sinusoids
            double lincom =0;
            /* Whole time series */
            for (i =0; i < NUM_POINTS; ++i) {
                /* Expansion term */
                double t = time_series[i];
                for (k =1; k < EXPAN_TERM; ++k) {
                    lincom += pow(-1.0,k+1) * sin(2*(k*t))/(k);
                }
                sawtooth_wave[i] = 0.5 + 2*lincom /(M_PI);
                // reset the accumulation term
                lincom =0;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    gnuplot_ctrl *h;
    /* Initialize the gnuplot handle */
    printf("*** DSP example of FFT gnuplot through C ***\n") ;
    h = gnuplot_init() ;

    if (argc >1) {
        EXPAN_TERM = atoi(argv[1]);
    } else if(argc ==3) {
        SLEEP_LGTH = atoi(argv[2]);
    }

    /* Initialize time series */
    double time_range[2] = {-4.0,4.0};
    init_time_series(time_range);

    // Fig.1
    fourier_expansion(SQUARE_WAVE);
    gnuplot_resetplot(h) ; gnuplot_setstyle(h, "lines") ;
    printf("\n\n*** Square Wave \n") ;
    gnuplot_plot_xy(h, time_series, square_wave, NUM_POINTS, "Square Wave") ;
    sleep(SLEEP_LGTH) ;

    // Fig.2
    fourier_expansion(TRIANGLE_WAVE);
    gnuplot_resetplot(h) ; gnuplot_setstyle(h, "lines") ;
    printf("\n\n*** Triangular Wave \n") ;
    gnuplot_plot_xy(h, time_series, triangle_wave, NUM_POINTS, "Triangular Wave") ;
    sleep(SLEEP_LGTH) ;

    // Fig.3
    fourier_expansion(SAWTOOTH_WAVE);
    gnuplot_resetplot(h) ; gnuplot_setstyle(h, "lines") ;
    printf("\n\n*** Sawtooth Wave \n") ;
    gnuplot_plot_xy(h, time_series, sawtooth_wave, NUM_POINTS, "Sawtooth Wave") ;
    sleep(SLEEP_LGTH) ;

    printf("\n\n*** end of DSP example ***\n") ;
    gnuplot_close(h);
    return 0;
}
