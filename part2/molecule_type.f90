module molecule_type
    use atom_type
    implicit none

    type molecule

      ! Private attributes
      type(atom), dimension(:), allocatable, public :: atoms
      integer , public :: nb_atoms
      logical, public :: with_bounds
      
      contains
        procedure :: init_mol, add_atom, print_mol

        ! allow calls like print *
        generic :: write(formatted) => print_mol

   end type molecule

  contains

    subroutine init_mol(m, size)

      ! Declaration
      class (molecule),intent(inout) :: m
      integer, intent(in) :: size

      ! Allocation
      allocate(m%atoms(size))

      ! Initialization
      m%nb_atoms=0

    end subroutine init_mol

    subroutine add_atom(m, a)

      ! Declaration
      class(molecule),intent(inout) :: m
      type(atom), intent(in) ::a 

      ! Enough place ?
      if(m%nb_atoms==size(m%atoms)) then
        ! No
        print *,"Not added atom: ", a
      else
        ! Yes ==> adding
        m%nb_atoms = m%nb_atoms+1
        m%atoms(m%nb_atoms)=a
      end if

    end subroutine add_atom

    subroutine print_mol(m, unit, iotype, v_list, iostat, iomsg)

      ! Declaration
      class (molecule),intent(in) :: m
      integer, intent(in) :: unit
      character(len=*), intent(in) :: iotype
      integer, dimension(:), intent(in):: v_list
      integer, intent(out) :: iostat
      character(len=*), intent(inout) :: iomsg

      ! Printing ==> call the home-made printing for atom_type    
      print *, m%atoms

      ! Everything ok here <=> need by ifort not gfortran
      iostat=0
    
    end subroutine print_mol

end module molecule_type
