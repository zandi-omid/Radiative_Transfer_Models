! Batch interface to the supplied fixed-form DMIESS routine.
program mie_batch
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: it = 183, ll = 9000
  real(8) :: theta(it), eltrmx(4,it,2), pi_angles(3,it), tau_angles(3,it)
  real(8) :: costheta(it), sin2theta(it), x, real_index, absorption_index
  real(8) :: wavelength, radius, wavenumber, qext, qscat, gqscat, qbs
  real(8) :: phase_forward, phase_backward, normalization
  complex(8) :: acap(ll)
  integer :: status, jx, j
  character(len=16) :: mode
  logical :: phase_mode

  call get_command_argument(1, mode)
  phase_mode = trim(mode) == '--phase'
  ! The default efficiency-only interface still uses theta=0.
  theta = 0.0d0
  do
    read(*, *, iostat=status) x, real_index, absorption_index, wavelength
    if (status < 0) exit
    if (status /= 0) stop 1
    if (x <= 0 .or. real_index <= 0 .or. absorption_index < 0 .or. wavelength <= 0) stop 2
    jx = 1
    if (phase_mode) then
      read(*, *, iostat=status) jx
      if (status /= 0 .or. jx < 1 .or. jx > it) stop 4
      read(*, *, iostat=status) theta(1:jx)
      if (status /= 0) stop 4
      if (any(theta(1:jx) < 0.0d0) .or. any(theta(1:jx) > 90.0d0)) stop 4
    end if
    wavenumber = 2.0d0 * acos(-1.0d0) / wavelength
    radius = x / wavenumber
    call dmiess(radius, real_index, absorption_index, theta, jx, &
                qext, qscat, gqscat, eltrmx, pi_angles, tau_angles, &
                costheta, sin2theta, acap, qbs, it, ll, &
                radius, real_index, absorption_index, wavenumber)
    if (.not. ieee_is_finite(qext) .or. .not. ieee_is_finite(qscat) .or. &
        .not. ieee_is_finite(gqscat) .or. qscat <= 0.0d0) stop 3
    write(*, '(7(ES25.16E3,1X))') x, radius, qext, qscat, qext-qscat, &
                                 gqscat/qscat, qscat/qext
    if (phase_mode) then
      ! Original ANRM = 2*pi/(wavenumber^2*Cscat) = 2/(x^2*Qscat).
      normalization = 2.0d0 / (x*x*qscat)
      do j = 1, jx
        phase_forward = normalization * (eltrmx(1,j,1) + eltrmx(2,j,1))
        phase_backward = normalization * (eltrmx(1,j,2) + eltrmx(2,j,2))
        if (.not. ieee_is_finite(phase_forward) .or. &
            .not. ieee_is_finite(phase_backward)) stop 5
        if (phase_forward < 0.0d0 .or. phase_backward < 0.0d0) stop 5
        write(*, '(3(ES25.16E3,1X))') theta(j), phase_forward, phase_backward
      end do
    end if
  end do
end program mie_batch
