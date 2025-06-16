module ao_integrals
  ! Implements one- and two-electron AO integrals via Obara-Saika recursion
  use basis_module
  use constants_module
  implicit none
  private
  public :: compute_overlap_OS, compute_eri_OS

contains

  subroutine compute_overlap_OS(shellA, shellB, S)
    ! Compute overlap integrals between two shells via Obara-Saika
    type(shell_t), intent(in) :: shellA, shellB
    real(dp), intent(out) :: S(:,:)
    integer :: iA, iB, x, y, z, nx, ny, nz
    ! Loop over primitives
    real(dp) :: prefac, P(3), AB2, gamma
    integer :: p, q
    ! Initialize S to zero
    S = 0.0_dp
    do p = 1, shellA%nprims
      do q = 1, shellB%nprims
        ! Compute combined exponent and center
        gamma = shellA%prims(p)%exponent + shellB%prims(q)%exponent
        P = (shellA%prims(p)%exponent*shellA%center + shellB%prims(q)%exponent*shellB%center) / gamma
        AB2 = sum((shellA%center - shellB%center)**2)
        prefac = shellA%prims(p)%coefficient * shellB%prims(q)%coefficient * &
                (pi/gamma)**1.5 * exp(-shellA%prims(p)%exponent*shellB%prims(q)%exponent/gamma*AB2)
        ! Recursion for each Cartesian component
        do x = 0, shellA%l + shellB%l
          do y = 0, shellA%l + shellB%l
            do z = 0, shellA%l + shellB%l
              S(x+1,y+1,z+1) = S(x+1,y+1,z+1) + prefac * OS_recursion(x, shellA%l, shellB%l, P(1)-shellA%center(1), P(1)-shellB%center(1), gamma)
              ! repeat for y and z components (not expanded here)
            end do
          end do
        end do
      end do
    end do
  end subroutine compute_overlap_OS

  pure function OS_recursion(n, la, lb, PA, PB, gamma) result(val)
    ! Obara-Saika recursion kernel for one dimension
    integer, intent(in) :: n, la, lb
    real(dp), intent(in) :: PA, PB, gamma
    real(dp) :: val
    ! Base case
    if (la == 0 .and. lb == 0) then
      if (n == 0) then
        val = 1.0_dp
      else
        val = 0.0_dp
      end if
    else
      ! Recursive relation (simplified placeholder)
      val = PA * OS_recursion(n-1, la-1, lb, PA, PB, gamma) + PB * OS_recursion(n-1, la, lb-1, PA, PB, gamma)
    end if
  end function OS_recursion

  subroutine compute_eri_OS(shellA, shellB, shellC, shellD, ERI)
    ! Compute two-electron integrals (ij|kl) via Obara-Saika
    type(shell_t), intent(in) :: shellA, shellB, shellC, shellD
    real(dp), intent(out) :: ERI(:,:,:,:)
    ! To implement: 4-index recursion over primitives and shells
    ! Follow standard OS formulae for 4-center integrals
  end subroutine compute_eri_OS

end module ao_integrals