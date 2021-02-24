/******************************************************************************
*  OTFFT Laptime Command by C99 Interface
******************************************************************************/
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "otfft_c.h"
#include "../stopwatch.h"

static inline int imin(const int x, const int y) { return x < y ? x : y; }

int main()
{
    const int n_max  = 22;
    const int N_max  = 1 << 22;
    const int nN_max = n_max*N_max;
    int i, n, N, LOOPS, p;
    double t;
    ccomplex_t* x;
    void* obj;
    counter_t t1, t2;

    setbuf(stdout, NULL);
    x = (ccomplex_t  *) simd_malloc(N_max*sizeof(ccomplex_t ));
    for (n = 1; n <= n_max; n++) {
        printf("2^(%2d):", n);
        N = 1 << n;
        LOOPS = imin(70, n*4) * (nN_max/(n*N));
        for (p = 0; p < N; p++) {
            t = (double)p / N;
            x[p].re = 10 * cos(3 * 2 * M_PI * t * t);
            x[p].im = 10 * sin(3 * 2 * M_PI * t * t);
        };
        obj = otfft_fft_new(N);
        t1 = get_counter();
        for (i = 0; i < LOOPS; i++) {
            otfft_fft_fwd(obj, x);
            otfft_fft_inv(obj, x);
        }
        t2 = get_counter();
        otfft_fft_delete(obj);
        printf("%11.2f[usec]\n", usec(t2 - t1) / LOOPS);
    }
    simd_free(x);
    return 0;
}
