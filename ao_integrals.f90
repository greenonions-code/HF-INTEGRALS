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

    ! Build recurrence for A-side
    do ai = 1, la
      do n = ai, maxL
        R(n) = PAx * R(n-1) + (n > 1) * ((n-1) / (2.0_dp * gamma)) * R(n-2)
      end do
    end do
    ! Build recurrence for B-side
    do bi = 1, lb
      do n = bi, maxL
        R(n) = PBx * R(n-1) + (n > 1) * ((n-1) / (2.0_dp * gamma)) * R(n-2)
      end do
    end do
  end subroutine compute_OS_1D

  subroutine compute_eri_OS(shellA, shellB, shellC, shellD, ERI)
    ! Placeholder for optimized 4-center integrals using similar iterative OS
    type(shell_t), intent(in) :: shellA, shellB, shellC, shellD
    real(dp), intent(out)     :: ERI(:,:,:,:)
    ! Implementation would mirror compute_overlap_OS with 3D->6D loops,
    ! screening, and reuse of 1D recurrences for all four shell pairs.
  end subroutine compute_eri_OS

end module ao_integrals