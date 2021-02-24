/******************************************************************************
*  OTFFT C99 Sample
******************************************************************************/

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include "otfft.h"

typedef double complex dcomplex;
static const dcomplex j = _Complex_I;

void rdft_fwd0(int N, const double *x, dcomplex *y)
{
    const double theta0 = 2*M_PI/N;
    for (int k = 0; k < N; k++) {
        dcomplex sum = 0;
        for (int p = 0; p < N; p++) {
            const double theta = k*p * theta0;
            const dcomplex wkp = cos(theta) - j*sin(theta);
            sum += x[p] * wkp;
        }
        y[k] = sum;
    }
}
void rdft_fwd(int N, const double *x, dcomplex *y)
{
    rdft_fwd0(N, x, y);
    for (int k = 0; k < N; k++) y[k] /= N;
}
void rdft_fwdn(int N, const double *x, dcomplex *y) { rdft_fwd(N, x, y); }

void check_rfft()
{
    static const int N = 256;
    double err;
    double *x = (double *) simd_malloc(N*sizeof(double));
    double *y = (double *) simd_malloc(N*sizeof(double));
    double *z = (double *) simd_malloc(N*sizeof(double));
    dcomplex *ox = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    dcomplex *oy = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    for (int p = 0; p < N; p++) z[p] = rand()%100;
    void *obj = otfft_rfft_new(N);
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    rdft_fwd(N, x, ox);
    /*************************************************************************/
    otfft_rfft_fwd(obj, y, oy);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = ox[k] - oy[k];
        err += cabs(d)*cabs(d);
    }
    printf("RFFT fwd: %2ld\n", lrint(log10(err)));
    otfft_rfft_inv(obj, oy, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("RFFT inv: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    rdft_fwd0(N, x, ox);
    /*************************************************************************/
    otfft_rfft_fwd0(obj, y, oy);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = ox[k] - oy[k];
        err += cabs(d)*cabs(d);
    }
    printf("RFFT fwd0: %2ld\n", lrint(log10(err)));
    otfft_rfft_invn(obj, oy, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("RFFT invn: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    rdft_fwdn(N, x, ox);
    /*************************************************************************/
    otfft_rfft_fwdn(obj, y, oy);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = ox[k] - oy[k];
        err += cabs(d)*cabs(d);
    }
    printf("RFFT fwdn: %2ld\n", lrint(log10(err)));
    otfft_rfft_inv0(obj, oy, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("RFFT inv0: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    otfft_rfft_delete(obj);
    simd_free(oy);
    simd_free(ox);
    simd_free(z);
    simd_free(y);
    simd_free(x);
}

int main()
{
    srand(time(NULL));
    check_rfft();
    return 0;
}
