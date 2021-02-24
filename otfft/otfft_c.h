/******************************************************************************
*  OTFFT C Header Version 11.5e
*
*  Copyright (c) 2016 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_c_h
#define otfft_c_h

#if __cplusplus
extern "C" {
#endif

#ifdef OTFFT_LIBRARY_BUILD
#include <otfft_exports.h>
#else
#define OTFFT_EXPORT
#endif

#include <math.h>

	typedef struct {
		double re, im;
	} ccomplex_t;

	//=============================================================================

	OTFFT_EXPORT void* simd_malloc(size_t n);
	OTFFT_EXPORT void simd_free(void* p);

	//=============================================================================

	OTFFT_EXPORT void* otfft_fft_new(int N);
	OTFFT_EXPORT void otfft_fft_delete(void* p);

	OTFFT_EXPORT void otfft_fft_fwd(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_fwd0(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_fwdu(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_fwdn(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_inv(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_inv0(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_invu(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_fft_invn(void* p, ccomplex_t* x);

	//=============================================================================

	OTFFT_EXPORT void* otfft_fft0_new(int N);
	OTFFT_EXPORT void otfft_fft0_delete(void* p);

	OTFFT_EXPORT void otfft_fft0_fwd(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_fwd0(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_fwdu(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_fwdn(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_inv(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_inv0(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_invu(void* p, ccomplex_t* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_fft0_invn(void* p, ccomplex_t* x, ccomplex_t* y);

	//=============================================================================

	OTFFT_EXPORT void* otfft_rfft_new(int N);
	OTFFT_EXPORT void otfft_rfft_delete(void* p);

	OTFFT_EXPORT void otfft_rfft_fwd(void* p, const double* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_rfft_fwd0(void* p, const double* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_rfft_fwdu(void* p, const double* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_rfft_fwdn(void* p, const double* x, ccomplex_t* y);
	OTFFT_EXPORT void otfft_rfft_inv(void* p, ccomplex_t* x, double* y);
	OTFFT_EXPORT void otfft_rfft_inv0(void* p, ccomplex_t* x, double* y);
	OTFFT_EXPORT void otfft_rfft_invu(void* p, ccomplex_t* x, double* y);
	OTFFT_EXPORT void otfft_rfft_invn(void* p, ccomplex_t* x, double* y);

	//=============================================================================

	OTFFT_EXPORT void* otfft_dct_new(int N);
	OTFFT_EXPORT void otfft_dct_delete(void* p);

	OTFFT_EXPORT void otfft_dct_fwd(void* p, double* x);
	OTFFT_EXPORT void otfft_dct_fwd0(void* p, double* x);
	OTFFT_EXPORT void otfft_dct_fwdn(void* p, double* x);
	OTFFT_EXPORT void otfft_dct_inv(void* p, double* x);
	OTFFT_EXPORT void otfft_dct_inv0(void* p, double* x);
	OTFFT_EXPORT void otfft_dct_invn(void* p, double* x);

	//=============================================================================

	OTFFT_EXPORT void* otfft_dct0_new(int N);
	OTFFT_EXPORT void otfft_dct0_delete(void* p);

	OTFFT_EXPORT void otfft_dct0_fwd(void* p, double* x, double* y, ccomplex_t* z);
	OTFFT_EXPORT void otfft_dct0_fwd0(void* p, double* x, double* y, ccomplex_t* z);
	OTFFT_EXPORT void otfft_dct0_fwdn(void* p, double* x, double* y, ccomplex_t* z);
	OTFFT_EXPORT void otfft_dct0_inv(void* p, double* x, double* y, ccomplex_t* z);
	OTFFT_EXPORT void otfft_dct0_inv0(void* p, double* x, double* y, ccomplex_t* z);
	OTFFT_EXPORT void otfft_dct0_invn(void* p, double* x, double* y, ccomplex_t* z);

	//=============================================================================

	OTFFT_EXPORT void* otfft_bluestein_new(int N);
	OTFFT_EXPORT void otfft_bluestein_delete(void* p);

	OTFFT_EXPORT void otfft_bluestein_fwd(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_fwd0(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_fwdu(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_fwdn(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_inv(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_inv0(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_invu(void* p, ccomplex_t* x);
	OTFFT_EXPORT void otfft_bluestein_invn(void* p, ccomplex_t* x);

	//=============================================================================

#if __cplusplus
}
#endif

#endif // otfft_c_h
