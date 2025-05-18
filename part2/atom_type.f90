module atom_type

   type atom

      ! Attributes
      integer :: number
      character(len=4), public :: element
      real, dimension(3), public :: coordinates
      character(len=6), public :: type
      real, public :: radius1, radius2, radius3
      integer, public :: subst_id
      character(len=6), public :: subst_name
      real, public :: charge
      character(len=10), public :: status_bit

   contains
      procedure :: print_atom, init_atom

      ! allow calls like print * instead of call print_atom
      generic :: write(formatted) => print_atom

   end type atom

contains

   ! Constructor
   subroutine init_atom(a, i, elt, coordinates, t)

      ! Declaration
      class (atom),intent(inout) :: a
      integer, intent(in) :: i
      character(len=4), intent(in)::elt
      real, dimension(3), intent(in) :: coordinates
      character(len=6), intent(in) :: t

      ! Initialization
      a%element=elt
      a%coordinates=coordinates
      a%type=t
      a%number=i

   end subroutine init_atom

   ! The signature has to compatible with the classic formatted print method
   subroutine print_atom(a, unit, iotype, v_list, iostat, iomsg)

      ! Declaration
      class (atom),intent(in) :: a

      integer, intent(in) :: unit
      character(len=*), intent(in) :: iotype
      integer, dimension(:), intent(in):: v_list
      integer, intent(out) :: iostat
      character(len=*), intent(inout) :: iomsg

      ! Generic printing
      ! the '/' character add a new line
      write(unit=unit, fmt='(i8,3x,a4,x,3f8.4,x,a6/)', iostat=iostat, iomsg=iomsg) a%number,a%element,a%coordinates,a%type

   end subroutine print_atom

end module atom_type
