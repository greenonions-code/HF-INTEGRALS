module ao_integrals
  ! Performance-optimized AO integrals via iterative Obara–Saika
  use basis_module
  use constants_module
  implicit none
  private
  public :: compute_overlap_OS, compute_eri_OS
  integer, parameter :: dp = kind(1.0d0)

contains

  subroutine compute_overlap_OS(shellA, shellB, S)
    ! Computes overlap integrals S_xyz for two shells using iterative OS
    type(shell_t), intent(in)    :: shellA, shellB
    real(dp), intent(out)        :: S(:,:,:)
    integer                       :: la, lb, p, q, nx, ny, nz, maxL
    real(dp)                      :: gamma, prefac, AB2
    real(dp)                      :: P(3), PA(3), PB(3)
    real(dp), allocatable         :: Rx(:), Ry(:), Rz(:)

    la = shellA%l; lb = shellB%l
    maxL = la + lb
    allocate(Rx(0:maxL), Ry(0:maxL), Rz(0:maxL))
    S = 0.0_dp

    do p = 1, shellA%nprims
      do q = 1, shellB%nprims
        gamma = shellA%prims(p)%exponent + shellB%prims(q)%exponent
        P = (shellA%prims(p)%exponent*shellA%center + shellB%prims(q)%exponent*shellB%center)/gamma
        PA = P - shellA%center; PB = P - shellB%center
        AB2 = sum((shellA%center - shellB%center)**2)
        prefac = shellA%prims(p)%coefficient * shellB%prims(q)%coefficient * &
                (pi/gamma)**1.5 * exp(-shellA%prims(p)%exponent*shellB%prims(q)%exponent/gamma*AB2)

        call compute_OS_1D(la, lb, PA(1), PB(1), gamma, Rx)
        call compute_OS_1D(la, lb, PA(2), PB(2), gamma, Ry)
        call compute_OS_1D(la, lb, PA(3), PB(3), gamma, Rz)

        do nx = 0, maxL
          do ny = 0, maxL
            do nz = 0, maxL
              S(nx+1,ny+1,nz+1) = S(nx+1,ny+1,nz+1) + prefac * Rx(nx) * Ry(ny) * Rz(nz)
            end do
          end do
        end do
      end do
    end do

    deallocate(Rx, Ry, Rz)
  end subroutine compute_overlap_OS

  subroutine compute_OS_1D(la, lb, PAx, PBx, gamma, R)
    ! Iterative one-dimensional OS recurrence (no function calls)
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

  subroutine compute_eri_OS(shellA, shellB, shellC, shellD, ERI)
    ! Compute two-electron integrals (ij|kl) via Obara–Saika recursion
    type(shell_t), intent(in) :: shellA, shellB, shellC, shellD
    real(dp), intent(out)     :: ERI(:,:,:,:)
    integer                    :: la, lb, lc, ld
    integer                    :: p, q, r, s
    integer                    :: nxAB, nyAB, nzAB, nxCD, nyCD, nzCD
    integer                    :: maxAB, maxCD, mmax
    real(dp)                   :: gammaAB, gammaCD, gamma, prefacAB, prefacCD, Tval
    real(dp)                   :: PA(3), PB(3), PC(3), PD(3), PAB(3), PCD(3)
    real(dp), allocatable      :: RABx(:), RABx2(:), RyAB(:), RzAB(:)
    real(dp), allocatable      :: RCDx(:), RyCD(:), RzCD(:)
    real(dp)                   :: BoysF

    la = shellA%l; lb = shellB%l
    lc = shellC%l; ld = shellD%l
    maxAB = la + lb; maxCD = lc + ld
    mmax = maxAB + maxCD

    allocate(RABx(0:maxAB), RyAB(0:maxAB), RzAB(0:maxAB))
    allocate(RCDx(0:maxCD), RyCD(0:maxCD), RzCD(0:maxCD))
    ERI = 0.0_dp

    do p = 1, shellA%nprims
      do q = 1, shellB%nprims
        gammaAB = shellA%prims(p)%exponent + shellB%prims(q)%exponent
        PAB = (shellA%prims(p)%exponent*shellA%center + shellB%prims(q)%exponent*shellB%center)/gammaAB
        PA = PAB - shellA%center; PB = PAB - shellB%center
        prefacAB = shellA%prims(p)%coefficient * shellB%prims(q)%coefficient * (pi/gammaAB)**1.5
        call compute_OS_1D(la, lb, PA(1), PB(1), gammaAB, RABx)
        call compute_OS_1D(la, lb, PA(2), PB(2), gammaAB, RyAB)
        call compute_OS_1D(la, lb, PA(3), PB(3), gammaAB, RzAB)

        do r = 1, shellC%nprims
          do s = 1, shellD%nprims
            gammaCD = shellC%prims(r)%exponent + shellD%prims(s)%exponent
            PCD = (shellC%prims(r)%exponent*shellC%center + shellD%prims(s)%exponent*shellD%center)/gammaCD
            PC = PCD - shellC%center; PD = PCD - shellD%center
            prefacCD = shellC%prims(r)%coefficient * shellD%prims(s)%coefficient * (pi/gammaCD)**1.5
            call compute_OS_1D(lc, ld, PC(1), PD(1), gammaCD, RCDx)
            call compute_OS_1D(lc, ld, PC(2), PD(2), gammaCD, RyCD)
            call compute_OS_1D(lc, ld, PC(3), PD(3), gammaCD, RzCD)

            gamma = gammaAB * gammaCD / (gammaAB + gammaCD)
            Tval = gamma * sum((PAB - PCD)**2)

            do m = 0, mmax
              BoysF = boys_function(m, Tval)
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

  pure function boys_function(m, T) result(Fm)
    ! Compute Boys function F_m(T)
    integer, intent(in) :: m
    real(dp), intent(in) :: T
    real(dp) :: Fm
    ! Placeholder: implement accurate Boys function via incomplete gamma or recurrence
    if (T < 1e-8_dp) then
      Fm = 1.0_dp/(2.0_dp*m+1.0_dp)
    else
      Fm = 0.5_dp * T**(-m-0.5_dp) * gamma(m+0.5_dp) * erf(sqrt(T))
    end if
  end function boys_function

end module ao_integrals
