module basis_module
  ! Defines data structures and routines for Gaussian basis sets
  implicit none
  private
  public :: primitive_t, shell_t, basis_set, load_basis, normalize_basis

  integer, parameter :: dp = kind(1.0d0)

  ! Primitive Gaussian type
  type :: primitive_t
    real(dp) :: exponent      ! Gaussian exponent (alpha)
    real(dp) :: coefficient   ! Contraction coefficient (d)
    real(dp) :: norm_const    ! Normalization constant
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
    integer :: nbf            ! total number of basis functions (cartesian)
  end type basis_set

contains

  subroutine load_basis(basis, filename)
    ! Load basis set parameters from a file and validate
    type(basis_set), intent(out) :: basis
    character(len=*), intent(in) :: filename
    integer :: i, j, ios, unit
    character(len=32) :: fmt

    ! Open file
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      stop 'Error opening basis file: '//filename
    end if

    ! Read number of shells
    read(unit, *, iostat=ios) basis%nshells
    if (ios /= 0 .or. basis%nshells <= 0) then
      stop 'Invalid number of shells'
    end if

    allocate(basis%shells(basis%nshells))
    basis%nbf = 0

    do i = 1, basis%nshells
      ! Read shell header: l, nprims, center
      read(unit, *, iostat=ios) basis%shells(i)%l, basis%shells(i)%nprims, basis%shells(i)%center
      if (ios /= 0 .or. basis%shells(i)%l < 0 .or. basis%shells(i)%nprims <= 0) then
        stop 'Invalid shell parameters at shell '//trim(adjustl(itoa(i)))
      end if

      allocate(basis%shells(i)%prims(basis%shells(i)%nprims))
      do j = 1, basis%shells(i)%nprims
        read(unit, *, iostat=ios) basis%shells(i)%prims(j)%exponent, basis%shells(i)%prims(j)%coefficient
        if (ios /= 0 .or. basis%shells(i)%prims(j)%exponent <= 0.0_dp) then
          stop 'Invalid primitive data at shell '//trim(adjustl(itoa(i)))//', prim '//trim(adjustl(itoa(j)))
        end if
      end do

      ! Update basis function count (cartesian functions)
      basis%nbf = basis%nbf + (basis%shells(i)%l + 1)*(basis%shells(i)%l + 2)/2
    end do
    close(unit)

    ! Normalize basis functions
    call normalize_basis(basis)
  end subroutine load_basis

  subroutine normalize_basis(basis)
    ! Compute normalization constants for primitives and contracted shells
    type(basis_set), intent(inout) :: basis
    integer :: i, j
    real(dp) :: lfac, prefac

    do i = 1, basis%nshells
      lfac = basis%shells(i)%l
      do j = 1, basis%shells(i)%nprims
        prefac = (2.0_dp* basis%shells(i)%prims(j)%exponent/pi)**(3.0_dp/2.0_dp)
        basis%shells(i)%prims(j)%norm_const = sqrt(prefac * (4.0_dp**lfac) * (basis%shells(i)%prims(j)%exponent**lfac) / &
             double(factorial2(2*lfac-1)) )
        ! Incorporate contraction coefficient
        basis%shells(i)%prims(j)%coefficient = basis%shells(i)%prims(j)%coefficient * basis%shells(i)%prims(j)%norm_const
      end do
    end do
  end subroutine normalize_basis

  !--------------------------!
  ! Utility functions below  !
  !--------------------------!

  recursive function factorial2(n) result(res)
    integer, intent(in) :: n
    integer :: res
    if (n <= 1) then
      res = 1
    else
      res = n * factorial2(n-2)
    end if
  end function factorial2

  pure function itoa(i) result(str)
    integer, intent(in) :: i
    character(len=12) :: str
    write(str, '(I0)') i
  end function itoa

end module basis_module