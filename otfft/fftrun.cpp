/******************************************************************************
*  FFT Tuning Command Version 11.5e
*
*  Copyright (c) 2021 S. Temerkhanov
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include "otfft_misc.h"
#include "otfft_avxdif4.h"
#include "otfft_avxdit4.h"
#include "otfft_avxdif8.h"
#include "otfft_avxdit8.h"
#include "otfft_avxdif16.h"
#include "otfft_avxdit16.h"
#include "otfft_sixstep.h"
#include "msleep.h"

#ifdef _MSC_VER
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;
using namespace OTFFT_MISC;

#define DELAY1 1
#define DELAY2 100
#define FACTOR 1

static const int n_max = 24;

template <class FFT>
double laptime(int loops, const FFT& fft, complex_t *x, complex_t *y)
{
    using namespace chrono;
    typedef microseconds::rep counter_t;
    counter_t ret = 0;
    const chrono::time_point t1 = high_resolution_clock::now();
    for (size_t j = 0; j < loops; j++)
    {
        fft.fwd(x, y);
        fft.inv(x, y);
    }
    const chrono::time_point t2 = high_resolution_clock::now();
    ret = (double)duration_cast<microseconds>(t2 - t1).count() / loops;

    return ret;
}

void initialize(const int N, complex_vector x)
{
    for (size_t p = 0; p < N; p++)
    {
        const double t = double(p)/N;
        x[p].Re = 10 * cos(3*2*M_PI*t*t);
        x[p].Im = 10 * sin(3*2*M_PI*t*t);
    }
}

double run(int type, int loops, int n)
{
    double lap;
    complex_vector x = (complex_vector) simd_malloc(n*sizeof(complex_t));
    complex_vector y = (complex_vector) simd_malloc(n*sizeof(complex_t));

    initialize(n, x);

    switch (type)
    {
    case 1:
        lap = laptime(loops, OTFFT_AVXDIF4::FFT0(n), x, y);
        break;
    case 2:
        lap = laptime(loops, OTFFT_AVXDIT4::FFT0(n), x, y);
        break;
    case 3:
        lap = laptime(loops, OTFFT_AVXDIF8::FFT0(n), x, y);
        break;
    case 4:
        lap = laptime(loops, OTFFT_AVXDIT8::FFT0(n), x, y);
        break;
    case 5:
        lap = laptime(loops, OTFFT_AVXDIF16::FFT0(n), x, y);
        break;
    case 6:
        lap = laptime(loops, OTFFT_AVXDIT16::FFT0(n), x, y);
        break;
    case 7:
        lap = laptime(loops, OTFFT_SixStep::FFT0(n), x, y);
        break;
    default:
        break;
    }

    simd_free(y);
    simd_free(x);

    return lap;
}

int main(int argc, char *argv[]) try
{
    if (argc == 4)
    {
        const int type = atoi(argv[1]);
        if (type <= 0 || type > 7)
            throw "Unknown type";
        const int loops = atoi(argv[2]);
        if (loops <= 0)
            throw "Invalid loop count";
        const int n = atoi(argv[3]);
        if (n <= 0)
            throw "Invalid FFT length";
        const double lap = run(type, loops, 1 << n);
        cout << "Lap time: " << scientific << lap << " us" << endl;
    }
    else if (argc == 3)
    {
        const int loops = atoi(argv[1]);
        if (loops <= 0)
            throw "Invalid loop count";
        const int n = atoi(argv[2]);
        if (n <= 0)
            throw "Invalid FFT length";
        for (int type = 1; type <= 7; type++)
        {
            double lap = run(type, loops, 1 << n);
            cout << type << ": Lap time: " << scientific << lap << " us" << endl;
        }
    }
    else throw "usage1: fftrun [type] [loops] [log2(length)]\n";
    return 0;
}
catch (const char* message) { cerr << message << endl; exit(1); }
catch (...) { cerr << "\n""*** exception! ***\n"; exit(1); }
