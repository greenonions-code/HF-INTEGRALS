module basis_module
  ! Defines data structures for Gaussian basis sets
  implicit none
  private
  public :: primitive_t, shell_t, basis_set

  integer, parameter :: dp = kind(1.0d0)

  ! Primitive Gaussian type
  type :: primitive_t
    real(dp) :: exponent      ! Gaussian exponent (alpha)
    real(dp) :: coefficient   ! Contraction coefficient (d)
  end type primitive_t

  ! Shell: contracted set of primitives
  type :: shell_t
    integer :: l              ! angular momentum quantum number (0 = s, 1 = p, etc.)
    integer :: nprims         ! number of primitives in this shell
    type(primitive_t), allocatable :: prims(:)  ! primitives array
    real(dp) :: center(3)     ! Cartesian coordinates of shell center
  end type shell_t

  ! Basis set collection of shells
  type :: basis_set
    integer :: nshells        ! total number of shells
    type(shell_t), allocatable :: shells(:)  ! array of shells
    integer :: nbf            ! total number of basis functions
  end type basis_set

contains

  subroutine load_basis(basis, filename)
    ! Load basis set parameters from a file
    type(basis_set), intent(out) :: basis
    character(len=*), intent(in) :: filename
    ! File format (example):
    ! nshells
    ! l  nprims  x  y  z
    ! alpha_1  coeff_1
    ! ...
    integer :: i, j, ios, unit
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       stop 'Error opening basis file'
    end if
    read(unit, *) basis%nshells
    allocate(basis%shells(basis%nshells))
    basis%nbf = 0
    do i = 1, basis%nshells
       read(unit, *) basis%shells(i)%l, basis%shells(i)%nprims, basis%shells(i)%center
       allocate(basis%shells(i)%prims(basis%shells(i)%nprims))
       do j = 1, basis%shells(i)%nprims
         read(unit, *) basis%shells(i)%prims(j)%exponent, basis%shells(i)%prims(j)%coefficient
       end do
       ! Update basis function count (cartesian functions)
       basis%nbf = basis%nbf + (basis%shells(i)%l + 1)*(basis%shells(i)%l + 2)/2
    end do
    close(unit)
  end subroutine load_basis

end module basis_module