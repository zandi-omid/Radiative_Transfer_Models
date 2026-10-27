! Batch interface to the supplied fixed-form DMIESS routine.
program mie_batch
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: it = 183, ll = 9000
  real(8) :: theta(it), eltrmx(4,it,2), pi_angles(3,it), tau_angles(3,it)
  real(8) :: costheta(it), sin2theta(it), x, real_index, absorption_index
  real(8) :: wavelength, radius, wavenumber, qext, qscat, gqscat, qbs
  complex(8) :: acap(ll)
  integer :: status

  ! Only theta=0 is needed here: integrated efficiencies/g are computed
  ! from the series coefficients, not numerical angle integration.
  theta = 0.0d0
  do
    read(*, *, iostat=status) x, real_index, absorption_index, wavelength
    if (status < 0) exit
    if (status /= 0) stop 1
    if (x <= 0 .or. real_index <= 0 .or. absorption_index < 0 .or. wavelength <= 0) stop 2
    wavenumber = 2.0d0 * acos(-1.0d0) / wavelength
    radius = x / wavenumber
    call dmiess(radius, real_index, absorption_index, theta, 1, &
                qext, qscat, gqscat, eltrmx, pi_angles, tau_angles, &
                costheta, sin2theta, acap, qbs, it, ll, &
                radius, real_index, absorption_index, wavenumber)
    if (.not. ieee_is_finite(qext) .or. .not. ieee_is_finite(qscat) .or. &
        .not. ieee_is_finite(gqscat) .or. qscat <= 0.0d0) stop 3
    write(*, '(7(ES25.16E3,1X))') x, radius, qext, qscat, qext-qscat, &
                                 gqscat/qscat, qscat/qext
  end do
end program mie_batch
