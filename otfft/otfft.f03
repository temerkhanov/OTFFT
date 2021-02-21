!==============================================================================
!  OTFFT Fortran Interface Version 11.5e
!
!  Copyright (c) 2016 OK Ojisan(Takuya OKAHISA)
!  Released under the MIT license
!  http://opensource.org/licenses/mit-license.php
!==============================================================================

interface
    !--------------------------------------------------------------------------

    function simd_malloc(n) bind(C, name='simd_malloc')
        use, intrinsic :: iso_c_binding
        integer(c_size_t), value :: n
        type(c_ptr) :: simd_malloc
    end function

    subroutine simd_free(p) bind(C, name='simd_free')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    !--------------------------------------------------------------------------

    function otfft_fft_new(n) bind(C, name='otfft_fft_new')
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n
        type(c_ptr) :: otfft_fft_new
    end function

    subroutine otfft_fft_delete(p) bind(C, name='otfft_fft_delete')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    subroutine otfft_fft_fwd(p, x) bind(C, name='otfft_fft_fwd')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_fft_fwd0(p, x) bind(C, name='otfft_fft_fwd0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_fft_fwdu(p, x) bind(C, name='otfft_fft_fwdu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_fft_fwdn(p, x) bind(C, name='otfft_fft_fwdn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine

    subroutine otfft_fft_inv(p, x) bind(C, name='otfft_fft_inv')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_fft_inv0(p, x) bind(C, name='otfft_fft_inv0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_fft_invu(p, x) bind(C, name='otfft_fft_invu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_fft_invn(p, x) bind(C, name='otfft_fft_invn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine

    !--------------------------------------------------------------------------

    function otfft_fft0_new(n) bind(C, name='otfft_fft0_new')
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n
        type(c_ptr) :: otfft_fft0_new
    end function

    subroutine otfft_fft0_delete(p) bind(C, name='otfft_fft0_delete')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    subroutine otfft_fft0_fwd(p, x, y) bind(C, name='otfft_fft0_fwd')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine
    subroutine otfft_fft0_fwd0(p, x, y) bind(C, name='otfft_fft0_fwd0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine
    subroutine otfft_fft0_fwdu(p, x, y) bind(C, name='otfft_fft0_fwdu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine
    subroutine otfft_fft0_fwdn(p, x, y) bind(C, name='otfft_fft0_fwdn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine

    subroutine otfft_fft0_inv(p, x, y) bind(C, name='otfft_fft0_inv')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine
    subroutine otfft_fft0_inv0(p, x, y) bind(C, name='otfft_fft0_inv0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine
    subroutine otfft_fft0_invu(p, x, y) bind(C, name='otfft_fft0_invu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine
    subroutine otfft_fft0_invn(p, x, y) bind(C, name='otfft_fft0_invn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        complex(c_double), dimension(*), intent(inout) :: y
    end subroutine

    !--------------------------------------------------------------------------

    function otfft_rfft_new(n) bind(C, name='otfft_rfft_new')
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n
        type(c_ptr) :: otfft_rfft_new
    end function

    subroutine otfft_rfft_delete(p) bind(C, name='otfft_rfft_delete')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    subroutine otfft_rfft_fwd(p, x, y) bind(C, name='otfft_rfft_fwd')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(in)  :: x
        complex(c_double), dimension(*), intent(out) :: y
    end subroutine
    subroutine otfft_rfft_fwd0(p, x, y) bind(C, name='otfft_rfft_fwd0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(in)  :: x
        complex(c_double), dimension(*), intent(out) :: y
    end subroutine
    subroutine otfft_rfft_fwdu(p, x, y) bind(C, name='otfft_rfft_fwdu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(in)  :: x
        complex(c_double), dimension(*), intent(out) :: y
    end subroutine
    subroutine otfft_rfft_fwdn(p, x, y) bind(C, name='otfft_rfft_fwdn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(in)  :: x
        complex(c_double), dimension(*), intent(out) :: y
    end subroutine

    subroutine otfft_rfft_inv(p, x, y) bind(C, name='otfft_rfft_inv')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(out)   :: y
    end subroutine
    subroutine otfft_rfft_inv0(p, x, y) bind(C, name='otfft_rfft_inv0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(out)   :: y
    end subroutine
    subroutine otfft_rfft_invu(p, x, y) bind(C, name='otfft_rfft_invu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(out)   :: y
    end subroutine
    subroutine otfft_rfft_invn(p, x, y) bind(C, name='otfft_rfft_invn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(out)   :: y
    end subroutine

    !--------------------------------------------------------------------------

    function otfft_dct_new(n) bind(C, name='otfft_dct_new')
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n
        type(c_ptr) :: otfft_dct_new
    end function

    subroutine otfft_dct_delete(p) bind(C, name='otfft_dct_delete')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    subroutine otfft_dct_fwd(p, x) bind(C, name='otfft_dct_fwd')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_dct_fwd0(p, x) bind(C, name='otfft_dct_fwd0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_dct_fwdn(p, x) bind(C, name='otfft_dct_fwdn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double), dimension(*), intent(inout) :: x
    end subroutine

    subroutine otfft_dct_inv(p, x) bind(C, name='otfft_dct_inv')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_dct_inv0(p, x) bind(C, name='otfft_dct_inv0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_dct_invn(p, x) bind(C, name='otfft_dct_invn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double), dimension(*), intent(inout) :: x
    end subroutine

    !--------------------------------------------------------------------------

    function otfft_dct0_new(n) bind(C, name='otfft_dct0_new')
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n
        type(c_ptr) :: otfft_dct0_new
    end function

    subroutine otfft_dct0_delete(p) bind(C, name='otfft_dct0_delete')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    subroutine otfft_dct0_fwd(p, x, y, z) bind(C, name='otfft_dct0_fwd')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(inout) :: y
        complex(c_double), dimension(*), intent(inout) :: z
    end subroutine
    subroutine otfft_dct0_fwd0(p, x, y, z) bind(C, name='otfft_dct0_fwd0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(inout) :: y
        complex(c_double), dimension(*), intent(inout) :: z
    end subroutine
    subroutine otfft_dct0_fwdn(p, x, y, z) bind(C, name='otfft_dct0_fwdn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(inout) :: y
        complex(c_double), dimension(*), intent(inout) :: z
    end subroutine

    subroutine otfft_dct0_inv(p, x, y, z) bind(C, name='otfft_dct0_inv')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(inout) :: y
        complex(c_double), dimension(*), intent(inout) :: z
    end subroutine
    subroutine otfft_dct0_inv0(p, x, y, z) bind(C, name='otfft_dct0_inv0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(inout) :: y
        complex(c_double), dimension(*), intent(inout) :: z
    end subroutine
    subroutine otfft_dct0_invn(p, x, y, z) bind(C, name='otfft_dct0_invn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        real(c_double),    dimension(*), intent(inout) :: x
        real(c_double),    dimension(*), intent(inout) :: y
        complex(c_double), dimension(*), intent(inout) :: z
    end subroutine

    !--------------------------------------------------------------------------

    function otfft_bluestein_new(n) bind(C, name='otfft_bluestein_new')
        use, intrinsic :: iso_c_binding
        integer(c_int), value :: n
        type(c_ptr) :: otfft_bluestein_new
    end function

    subroutine otfft_bluestein_delete(p) bind(C, name='otfft_bluestein_delete')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
    end subroutine

    subroutine otfft_bluestein_fwd(p, x) bind(C, name='otfft_bluestein_fwd')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_bluestein_fwd0(p, x) bind(C, name='otfft_bluestein_fwd0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_bluestein_fwdu(p, x) bind(C, name='otfft_bluestein_fwdu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_bluestein_fwdn(p, x) bind(C, name='otfft_bluestein_fwdn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine

    subroutine otfft_bluestein_inv(p, x) bind(C, name='otfft_bluestein_inv')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_bluestein_inv0(p, x) bind(C, name='otfft_bluestein_inv0')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_bluestein_invu(p, x) bind(C, name='otfft_bluestein_invu')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine
    subroutine otfft_bluestein_invn(p, x) bind(C, name='otfft_bluestein_invn')
        use, intrinsic :: iso_c_binding
        type(c_ptr), value :: p
        complex(c_double), dimension(*), intent(inout) :: x
    end subroutine

    !--------------------------------------------------------------------------
end interface

!==============================================================================
