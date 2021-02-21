/******************************************************************************
*  OTFFT C99 Sample
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include "otfft/otfft.h"

typedef double complex dcomplex;

void dct_fwd0(int N, double *x)
{
    const double theta0 = M_PI/N;
    double *y = (double *) simd_malloc(N*sizeof(double));
    for (int k = 0; k < N; k++) {
        double sum = 0;
        for (int p = 0; p < N; p++) {
            sum += x[p] * cos(theta0*k*(p+1.0/2));
        }
        y[k] = sum;
    }
    for (int k = 0; k < N; k++) x[k] = y[k];
    simd_free(y);
}
void dct_fwd(int N, double *x)
{
    dct_fwd0(N, x);
    for (int k = 0; k < N; k++) x[k] /= N;
}
void dct_fwdn(int N, double *x) { dct_fwd(N, x); }

void check_dct()
{
    static const int N = 256;
    double err;
    double *x = (double *) simd_malloc(N*sizeof(double));
    double *y = (double *) simd_malloc(N*sizeof(double));
    double *z = (double *) simd_malloc(N*sizeof(double));
    for (int p = 0; p < N; p++) z[p] = rand()%100;
    void *obj = otfft_dct_new(N);
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dct_fwd(N, x);
    /*************************************************************************/
    otfft_dct_fwd(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = x[k] - y[k];
        err += d*d;
    }
    printf("DCT fwd: %2ld\n", lrint(log10(err)));
    otfft_dct_inv(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("DCT inv: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dct_fwd0(N, x);
    /*************************************************************************/
    otfft_dct_fwd0(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = x[k] - y[k];
        err += d*d;
    }
    printf("DCT fwd0: %2ld\n", lrint(log10(err)));
    otfft_dct_invn(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("DCT invn: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dct_fwdn(N, x);
    /*************************************************************************/
    otfft_dct_fwdn(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = x[k] - y[k];
        err += d*d;
    }
    printf("DCT fwdn: %2ld\n", lrint(log10(err)));
    otfft_dct_inv0(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("DCT inv0: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    otfft_dct_delete(obj);
    simd_free(z);
    simd_free(y);
    simd_free(x);
}

void check_dct0()
{
    static const int N = 256;
    double err;
    double *x = (double *) simd_malloc(N*sizeof(double));
    double *y = (double *) simd_malloc(N*sizeof(double));
    double *z = (double *) simd_malloc(N*sizeof(double));
    double *v = (double *) simd_malloc(N*sizeof(double));
    dcomplex *w = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    for (int p = 0; p < N; p++) z[p] = rand()%100;
    void *obj = otfft_dct0_new(N);
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dct_fwd(N, x);
    /*************************************************************************/
    otfft_dct0_fwd(obj, y, v, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = x[k] - y[k];
        err += d*d;
    }
    printf("DCT0 fwd: %2ld\n", lrint(log10(err)));
    otfft_dct0_inv(obj, y, v, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("DCT0 inv: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dct_fwd0(N, x);
    /*************************************************************************/
    otfft_dct0_fwd0(obj, y, v, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = x[k] - y[k];
        err += d*d;
    }
    printf("DCT0 fwd0: %2ld\n", lrint(log10(err)));
    otfft_dct0_invn(obj, y, v, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("DCT0 invn: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dct_fwdn(N, x);
    /*************************************************************************/
    otfft_dct0_fwdn(obj, y, v, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = x[k] - y[k];
        err += d*d;
    }
    printf("DCT0 fwdn: %2ld\n", lrint(log10(err)));
    otfft_dct0_inv0(obj, y, v, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const double d = y[k] - z[k];
        err += d*d;
    }
    printf("DCT0 inv0: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    otfft_dct0_delete(obj);
    simd_free(w);
    simd_free(v);
    simd_free(z);
    simd_free(y);
    simd_free(x);
}

int main()
{
    srand(time(NULL));
    check_dct();
    printf("\n");
    check_dct0();
    return 0;
}
