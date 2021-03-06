!==============================================================================
!  OTFFT RFFT Laptime Command by Fortran Interface
!==============================================================================

program main
    use, intrinsic :: iso_c_binding
    implicit none
    include "otfft/otfft.f03"
    integer, parameter :: lnmax = 22, nmax = 2**lnmax
    integer, parameter :: sizeof_real = 8, sizeof_complex = 16
    double precision, parameter :: pi = 3.14159265358979323846264338327950288d0
    real(c_double), pointer :: x(:)
    complex(c_double), pointer :: y(:)
    !real(c_double), allocatable :: x(:)
    !complex(c_double), allocatable :: y(:)
    integer n, ln, b, p, LOOPS
    integer t1, t2, cps, tmax
    double precision t
    type(c_ptr) mem1, mem2, obj

    mem1 = simd_malloc(int(nmax*sizeof_real, c_size_t))
    mem2 = simd_malloc(int(nmax*sizeof_complex, c_size_t))
    call c_f_pointer(mem1, x, [nmax])
    call c_f_pointer(mem2, y, [nmax])
    !allocate(x(0:nmax-1), y(0:nmax-1))
    do ln = 1, lnmax
        print "('2^(',i2,'):',$)", ln
        n = 2**ln
        b = lbound(x, 1)
        do p = b, n-(1-b)
            t = dble(p - b) / n
            x(p) = 10 * real(cos(3*2*pi*t*t), c_double)
        end do
        LOOPS = min(70, n*4) * ((nmax*lnmax) / (n*ln))
        obj = otfft_rfft_new(int(n, c_int))
        call system_clock(t1, cps, tmax)
        do p = 1, LOOPS
            call otfft_rfft_fwd(obj, x, y)
            call otfft_rfft_inv(obj, y, x)
        end do
        call system_clock(t2)
        call otfft_rfft_delete(obj)
        print "(f11.2,'[usec]')", usec(t2 - t1) / LOOPS
    end do
    !deallocate(y, x)
    call simd_free(mem2)
    call simd_free(mem1)

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
