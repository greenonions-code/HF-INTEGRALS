module ao_integrals
  ! Performance-optimized AO integrals via iterative Obara–Saika with Boys caching
  use basis_module
  use constants_module
  implicit none
  private
  public :: init_boys_cache, compute_overlap_OS, compute_eri_OS, get_boys

  integer, parameter :: dp = kind(1.0d0)

  ! Boys caching parameters
  integer :: boys_mmax = -1
  integer :: boys_nT = -1
  real(dp) :: boys_Tmin, boys_Tmax, boys_dT
  real(dp), allocatable :: boys_table(:,:)

contains

  subroutine init_boys_cache(mmax, nT, Tmin, Tmax)
    ! Initialize Boys function cache for 0<=m<=mmax over nT grid points
    integer, intent(in) :: mmax, nT
    real(dp), intent(in) :: Tmin, Tmax
    integer :: m, i
    real(dp) :: T

    if (Tmax <= Tmin .or. nT < 2 .or. mmax < 0) then
      stop 'Invalid Boys cache parameters'
    end if
    boys_mmax = mmax
    boys_nT = nT
    boys_Tmin = Tmin
    boys_Tmax = Tmax
    boys_dT = (Tmax - Tmin) / real(nT-1, dp)
    allocate(boys_table(0:boys_mmax, 0:boys_nT-1))

    do m = 0, boys_mmax
      do i = 0, boys_nT-1
        T = boys_Tmin + real(i, dp) * boys_dT
        boys_table(m,i) = compute_boys_direct(m, T)
      end do
    end do
  end subroutine init_boys_cache

  function get_boys(m, T) result(Fm)
    ! Retrieve or interpolate Boys function value from cache, fallback if needed
    integer, intent(in) :: m
    real(dp), intent(in) :: T
    real(dp) :: Fm
    integer :: idx
    real(dp) :: tloc, w

    if (allocated(boys_table) .and. m >= 0 .and. m <= boys_mmax) then
      if (T <= boys_Tmin) then
        Fm = boys_table(m,0)
      else if (T >= boys_Tmax) then
        Fm = compute_boys_direct(m, T)
      else
        tloc = (T - boys_Tmin) / boys_dT
        idx = int(floor(tloc))
        w = tloc - real(idx, dp)
        Fm = (1.0_dp - w)*boys_table(m, idx) + w*boys_table(m, idx+1)
      end if
    else
      Fm = compute_boys_direct(m, T)
    end if
  end function get_boys

  pure function compute_boys_direct(m, T) result(Fm)
    ! Direct evaluation of Boys function F_m(T)
    integer, intent(in) :: m
    real(dp), intent(in) :: T
    real(dp) :: Fm
    if (T < 1e-8_dp) then
      Fm = 1.0_dp / (2.0_dp*real(m,dp) + 1.0_dp)
    else
      Fm = 0.5_dp * T**(-m-0.5_dp) * gamma(m+0.5_dp) * erf(sqrt(T))
    end if
  end function compute_boys_direct

  subroutine compute_overlap_OS(shellA, shellB, S)
    ! (unchanged) ...
    ! Implementation same as before
  end subroutine compute_overlap_OS

  subroutine compute_eri_OS(shellA, shellB, shellC, shellD, ERI)
    ! Compute two-electron integrals (ij|kl) via OS recursion using cached Boys
    type(shell_t), intent(in) :: shellA, shellB, shellC, shellD
    real(dp), intent(out)     :: ERI(:,:,:,:)
    integer                    :: la, lb, lc, ld, p, q, r, s
    integer                    :: nxAB, nyAB, nzAB, nxCD, nyCD, nzCD, maxAB, maxCD, mmax
    real(dp)                   :: gammaAB, gammaCD, gamma, prefacAB, prefacCD, Tval
    real(dp)                   :: PA(3), PB(3), PC(3), PD(3), PAB(3), PCD(3)
    real(dp), allocatable      :: RABx(:), RyAB(:), RzAB(:)
    real(dp), allocatable      :: RCDx(:), RyCD(:), RzCD(:)
    real(dp)                   :: BoysF

    la = shellA%l; lb = shellB%l
    lc = shellC%l; ld = shellD%l
    maxAB = la + lb; maxCD = lc + ld
    mmax = maxAB + maxCD

    ! Ensure cache covers required m and a reasonable T-range
    if (.not. allocated(boys_table) .or. mmax > boys_mmax) then
      call init_boys_cache(mmax, 100, 0.0_dp, 50.0_dp)
    end if

    allocate(RABx(0:maxAB), RyAB(0:maxAB), RzAB(0:maxAB))
    allocate(RCDx(0:maxCD), RyCD(0:maxCD), RzCD(0:maxCD))
    ERI = 0.0_dp

    do p = 1, shellA%nprims
      do q = 1, shellB%nprims
        gammaAB = shellA%prims(p)%exponent + shellB%prims(q)%exponent
        PAB = (shellA%prims(p)%exponent*shellA%center + shellB%prims(q)%exponent*shellB%center) / gammaAB
        PA = PAB - shellA%center; PB = PAB - shellB%center
        prefacAB = shellA%prims(p)%coefficient * shellB%prims(q)%coefficient * (pi/gammaAB)**1.5
        call compute_OS_1D(la, lb, PA(1), PB(1), gammaAB, RABx)
        call compute_OS_1D(la, lb, PA(2), PB(2), gammaAB, RyAB)
        call compute_OS_1D(la, lb, PA(3), PB(3), gammaAB, RzAB)

        do r = 1, shellC%nprims
          do s = 1, shellD%nprims
            gammaCD = shellC%prims(r)%exponent + shellD%prims(s)%exponent
            PCD = (shellC%prims(r)%exponent*shellC%center + shellD%prims(s)%exponent*shellD%center) / gammaCD
            PC = PCD - shellC%center; PD = PCD - shellD%center
            prefacCD = shellC%prims(r)%coefficient * shellD%prims(s)%coefficient * (pi/gammaCD)**1.5
            call compute_OS_1D(lc, ld, PC(1), PD(1), gammaCD, RCDx)
            call compute_OS_1D(lc, ld, PC(2), PD(2), gammaCD, RyCD)
            call compute_OS_1D(lc, ld, PC(3), PD(3), gammaCD, RzCD)

            gamma = gammaAB * gammaCD / (gammaAB + gammaCD)
            Tval = gamma * sum((PAB - PCD)**2)

            do m = 0, mmax
              BoysF = get_boys(m, Tval)
              do nxAB = 0, maxAB
                do nyAB = 0, maxAB
                  do nzAB = 0, maxAB
                    do nxCD = 0, maxCD
                      do nyCD = 0, maxCD
                        do nzCD = 0, maxCD
                          ERI(nxAB+1,nyAB+1,nzAB+1,nxCD+1,nyCD+1,nzCD+1) = &
                            ERI(nxAB+1,nyAB+1,nzAB+1,nxCD+1,nyCD+1,nzCD+1) + &
                            prefacAB*prefacCD * RABx(nxAB)*RyAB(nyAB)*RzAB(nzAB) * &
                            RCDx(nxCD)*RyCD(nyCD)*RzCD(nzCD) * BoysF
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    deallocate(RABx, RyAB, RzAB, RCDx, RyCD, RzCD)
  end subroutine compute_eri_OS

  ! compute_OS_1D remains unchanged
  subroutine compute_OS_1D(la, lb, PAx, PBx, gamma, R)
    integer, intent(in)    :: la, lb
    real(dp), intent(in)   :: PAx, PBx, gamma
    real(dp), intent(out)  :: R(0:)
    integer                 :: ai, bi, n, maxL
    maxL = la + lb
    R = 0.0_dp; R(0) = 1.0_dp
    do ai = 1, la
      do n = ai, maxL
        R(n) = PAx * R(n-1) + (n-1)/(2.0_dp*gamma) * ((n-1) * R(n-2))
      end do
    end do
    do bi = 1, lb
      do n = bi, maxL
        R(n) = PBx * R(n-1) + (n-1)/(2.0_dp*gamma) * ((n-1) * R(n-2))
      end do
    end do
  end subroutine compute_OS_1D

end module ao_integrals
