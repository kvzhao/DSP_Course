
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
double reconstructed[NUM_POINTS];

double time_series     [NUM_POINTS];
double square   [NUM_POINTS];
double skew     [NUM_POINTS];
double triangle [NUM_POINTS];

void init_square_wave(double t_range[2], double period) {
    double time_length = t_range[1] - t_range[0];
    double time_step = time_length / NUM_POINTS;
    double t;
    int i, A = 1;
    int pn = period / time_step, p2 = pn >> 2;

    for (i =0, t = t_range[0]; i< NUM_POINTS; i++, t += time_step) {
        if ((i%pn) > p2)
            square[i] = A;
        else
            square[i] = -A;
        time_series[i] = t;
    }
}

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

    }
}

/* Resume reading here */

int main() {

    gnuplot_ctrl *h;

    /* Initialize the gnuplot handle */
    printf("*** DSP example of FFT gnuplot through C ***\n") ;
    h = gnuplot_init() ;

/*
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
    */

    fftw_complex sqrs[NUM_POINTS/2 +1];

    double time_range[2] = {0.0,4.0};
    double time_period = 2;
    init_square_wave(time_range, time_period);


    unsigned int n = NUM_POINTS/ 2;
    fftw_plan plan;

    plan = fftw_plan_dft_r2c_1d(NUM_POINTS, square, sqrs, FFTW_ESTIMATE);
    fftw_execute(plan);

    plan = fftw_plan_dft_c2r_1d(NUM_POINTS, sqrs, reconstructed, FFTW_ESTIMATE);
    fftw_execute(plan);

    /* Plot the signal */
    gnuplot_resetplot(h) ;
    gnuplot_setstyle(h, "lines") ;
    printf("\n\n*** Original Signal \n") ;
    gnuplot_plot_xy(h, time_series, reconstructed, NUM_POINTS, "Signal") ;
    sleep(SLEEP_LGTH) ;

    printf("\n\n") ;
    printf("*** end of DSP example ***\n") ;
    gnuplot_resetplot(h);
    gnuplot_close(h);

    fftw_destroy_plan(plan);
    return 0;
}
