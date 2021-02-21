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
static const dcomplex j = _Complex_I;

void dft_fwd0(int N, dcomplex *x)
{
    const double theta0 = 2*M_PI/N;
    dcomplex *y = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    for (int k = 0; k < N; k++) {
        dcomplex sum = 0;
        for (int p = 0; p < N; p++) {
            const double theta = k*p * theta0;
            const dcomplex wkp = cos(theta) - j*sin(theta);
            sum += x[p] * wkp;
        }
        y[k] = sum;
    }
    for (int k = 0; k < N; k++) x[k] = y[k];
    simd_free(y);
}
void dft_fwd(int N, dcomplex *x)
{
    dft_fwd0(N, x);
    for (int k = 0; k < N; k++) x[k] /= N;
}
void dft_fwdn(int N, dcomplex *x) { dft_fwd(N, x); }

void check_fft()
{
    static const int N = 256;
    double err;
    dcomplex *x = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    dcomplex *y = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    dcomplex *z = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    for (int p = 0; p < N; p++) z[p] = rand()%100 + j*(rand()%100);
    void *obj = otfft_fft_new(N);
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dft_fwd(N, x);
    /*************************************************************************/
    otfft_fft_fwd(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = x[k] - y[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT fwd: %2ld\n", lrint(log10(err)));
    otfft_fft_inv(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = y[k] - z[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT inv: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dft_fwd0(N, x);
    /*************************************************************************/
    otfft_fft_fwd0(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = x[k] - y[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT fwd0: %2ld\n", lrint(log10(err)));
    otfft_fft_invn(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = y[k] - z[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT invn: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dft_fwdn(N, x);
    /*************************************************************************/
    otfft_fft_fwdn(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = x[k] - y[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT fwdn: %2ld\n", lrint(log10(err)));
    otfft_fft_inv0(obj, y);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = y[k] - z[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT inv0: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    otfft_fft_delete(obj);
    simd_free(z);
    simd_free(y);
    simd_free(x);
}

void check_fft0()
{
    static const int N = 256;
    double err;
    dcomplex *x = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    dcomplex *y = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    dcomplex *z = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    dcomplex *w = (dcomplex *) simd_malloc(N*sizeof(dcomplex));
    for (int p = 0; p < N; p++) z[p] = rand()%100 + j*(rand()%100);
    void *obj = otfft_fft0_new(N);
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dft_fwd(N, x);
    /*************************************************************************/
    otfft_fft0_fwd(obj, y, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = x[k] - y[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT0 fwd: %2ld\n", lrint(log10(err)));
    otfft_fft0_inv(obj, y, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = y[k] - z[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT0 inv: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dft_fwd0(N, x);
    /*************************************************************************/
    otfft_fft0_fwd0(obj, y, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = x[k] - y[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT0 fwd0: %2ld\n", lrint(log10(err)));
    otfft_fft0_invn(obj, y, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = y[k] - z[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT0 invn: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    for (int p = 0; p < N; p++) x[p] = y[p] = z[p];
    dft_fwdn(N, x);
    /*************************************************************************/
    otfft_fft0_fwdn(obj, y, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = x[k] - y[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT0 fwdn: %2ld\n", lrint(log10(err)));
    otfft_fft0_inv0(obj, y, w);
    err = 0;
    for (int k = 0; k < N; k++) {
        const dcomplex d = y[k] - z[k];
        err += cabs(d)*cabs(d);
    }
    printf("FFT0 inv0: %2ld\n", lrint(log10(err)));
    /*************************************************************************/
    otfft_fft0_delete(obj);
    simd_free(w);
    simd_free(z);
    simd_free(y);
    simd_free(x);
}

int main()
{
    srand(time(NULL));
    check_fft();
    printf("\n");
    check_fft0();
    return 0;
}
