/******************************************************************************
*  OTFFT Laptime Command by C99 Interface
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "otfft/otfft.h"
#include "../stopwatch.h"

typedef double _Complex dcomplex;
static const dcomplex j = _Complex_I;

inline int imin(const int x, const int y) { return x < y ? x : y; }

int main()
{
    static const int n_max  = 22;
    static const int N_max  = 1 << n_max;
    static const int nN_max = n_max*N_max;

    setbuf(stdout, NULL);
    dcomplex *x = (dcomplex *) simd_malloc(N_max*sizeof(dcomplex));
    for (int n = 1; n <= n_max; n++) {
        printf("2^(%2d):", n);
        const int N = 1 << n;
        const int LOOPS = imin(70, n*4) * (nN_max/(n*N));
        for (int p = 0; p < N; p++) {
            const double t = (double) p / N;
            x[p] = 10*cos(3*2*M_PI*t*t) + 10*j*sin(3*2*M_PI*t*t);
        }
        void *obj = otfft_fft_new(N);
        const counter_t t1 = get_counter();
        for (int i = 0; i < LOOPS; i++) {
            otfft_fft_fwd(obj, x);
            otfft_fft_inv(obj, x);
        }
        const counter_t t2 = get_counter();
        otfft_fft_delete(obj);
        printf("%11.2f[usec]\n", usec(t2 - t1) / LOOPS);
    }
    simd_free(x);
    return 0;
}
