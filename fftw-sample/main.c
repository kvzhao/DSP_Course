
/* Start reading here */

#include <fftw3.h>
#include "gnuplot_i.h"
#define NUM_POINTS 256


/* Never mind this bit */

#include <stdio.h>
#include <math.h>

#define REAL 0
#define IMAG 1
#define SLEEP_LGTH  2

double mag[NUM_POINTS];
double real_sig[NUM_POINTS];

void acquire_from_somewhere(fftw_complex* signal) {
    /* Generate two sine waves of different frequencies and
     * amplitudes.
     */

    int i;
    for (i = 0; i < NUM_POINTS; ++i) {
        double theta = (double)i / (double)NUM_POINTS * M_PI;

        signal[i][REAL] = 1.0 * cos(25.0 * theta) +
                          0.5 * cos(125.0 * theta);

        signal[i][IMAG] = 1.0 * sin(25.0 * theta) +
                          0.5 * sin(125.0 * theta);

        real_sig[i]     = signal[i][REAL];
    }
}

void do_something_with(fftw_complex* result) {
    int i;
    for (i = 0; i < NUM_POINTS; ++i) {
        mag[i] = sqrt(result[i][REAL] * result[i][REAL] +
                          result[i][IMAG] * result[i][IMAG]);
    }
}


/* Resume reading here */

int main() {

    gnuplot_ctrl *h;

    /* Initialize the gnuplot handle */
    printf("*** DSP example of FFT gnuplot through C ***\n") ;
    h = gnuplot_init() ;


    fftw_complex signal[NUM_POINTS];
    fftw_complex result[NUM_POINTS];

    fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS,
                                      signal,
                                      result,
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE);

    acquire_from_somewhere(signal);
    fftw_execute(plan);
    do_something_with(result);

    /* Plot the signal */
    gnuplot_resetplot(h) ;
    gnuplot_setstyle(h, "lines") ;
    printf("\n\n*** Original Signal \n") ;
    gnuplot_plot_x(h, real_sig, NUM_POINTS, "Signal") ;
    sleep(SLEEP_LGTH) ;
    gnuplot_resetplot(h);

    /* Plot the signal */
    gnuplot_resetplot(h) ;
    gnuplot_setstyle(h, "points") ;
    printf("\n\n*** FFTW Spectrum ***\n") ;
    gnuplot_plot_x(h, mag, NUM_POINTS, "FFT") ;
    sleep(SLEEP_LGTH) ;

    gnuplot_resetplot(h);
    printf("\n\n") ;
    printf("*** end of DSP example ***\n") ;
    gnuplot_close(h);
    fftw_destroy_plan(plan);

    return 0;
}
