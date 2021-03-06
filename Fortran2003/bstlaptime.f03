!==============================================================================
!  OTFFT Bluestein's FFT Laptime Command by Fortran Interface
!==============================================================================

program main
    use, intrinsic :: iso_c_binding
    implicit none
    include "otfft/otfft.f03"
    integer, parameter :: lnmax = 6, nmax = 10**lnmax
    integer, parameter :: sizeof_complex = 16
    double precision, parameter :: pi = 3.14159265358979323846264338327950288d0
    complex(c_double), pointer :: x(:)
    !complex(c_double), allocatable :: x(:)
    integer n, ln, b, p, LOOPS
    integer t1, t2, cps, tmax
    double precision t
    type(c_ptr) mem, obj

    mem = simd_malloc(int(nmax*sizeof_complex, c_size_t))
    call c_f_pointer(mem, x, [nmax])
    !allocate(x(0:nmax-1))
    do ln = 1, lnmax
        print "('10^(',i1,'):',$)", ln
        n = 10**ln
        b = lbound(x, 1)
        do p = b, n-(1-b)
            t = dble(p - b) / n
            x(p) = 10 * cmplx(cos(3*2*pi*t*t), sin(3*2*pi*t*t), c_double)
        end do
        LOOPS = min(10, n*3) * ((nmax*lnmax) / (n*ln))
        obj = otfft_bluestein_new(int(n, c_int))
        call system_clock(t1, cps, tmax)
        do p = 1, LOOPS
            call otfft_bluestein_fwd(obj, x)
            call otfft_bluestein_inv(obj, x)
        end do
        call system_clock(t2)
        call otfft_bluestein_delete(obj)
        print "(f11.2,'[usec]')", usec(t2 - t1) / LOOPS
    end do
    !deallocate(x);
    call simd_free(mem)

contains
    function usec(dt)
        integer, intent(in) :: dt
        double precision :: usec
        integer d

        if (dt >= 0) then
            d = dt
        else
            d = dt + tmax + 1
        end if
        usec = 1000000 * dble(d) / cps
    end function
end program
