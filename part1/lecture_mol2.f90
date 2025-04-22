program lecture_mol2
   use atom_type

   implicit none

   integer :: i, nbArg, end, ok
   character(len=128) :: fileName, line
   integer :: atomNumber, subStructureId
   character(len=4) :: atomName
   character(len=6) :: atomType
   character(len=6) :: subStructureName
   real :: x, y, z, charge
   real, dimension(3) :: coordinates
   character(len=80) :: mol_format = '(i5, a4, 3f8.4, a6, i5, a6, f8.4)'
   character(len=40) :: print_format = '(a4, 3f8.4)'

   integer :: unit = 10

   type(atom) :: a
   type(atom), dimension(:), allocatable :: atoms

   allocate(atoms(0))

   if(iargc() /= 1) then
      print '(a)', "Please provide a file name"
      stop unit
   end if

   call getarg(1, fileName)

   print '(a,a)', "File to read = ", trim(fileName)

   open(unit=unit, file=fileName, iostat=ok, status='old')

   if(ok /= 0) then
      print '(a,4x,a)', "Error during opening", fileName
      stop 20
   end if

   do
      read(unit, '(a)', iostat=end) line
      if (end /= 0) exit
      if (trim(line) == '@<TRIPOS>ATOM') then
         exit
      end if
   end do

   do
      !read(unit, '(a)', iostat=end) line
      read(unit, *, iostat=end) atomNumber, atomName, x, y, z, atomType, subStructureId, subStructureName, charge
      if (end /= 0) exit

      coordinates = (/x, y, z/)

      call a%init_atom(atomNumber, atomName, coordinates, atomType)

      call rayon_covalence(a)

      atoms = [atoms, a]
   end do

   close(unit)

   do i = 1, size(atoms)
      print *, atoms(i)
      print *, "Element:", trim(atoms(i)%element(1:1))
      print *, "Single Bond Radius:", atoms(i)%radius1
      print *, "Double Bond Radius:", atoms(i)%radius2
      print *, "Triple Bond Radius:", atoms(i)%radius3
   end do

   deallocate(atoms)

end program lecture_mol2

subroutine rayon_covalence(current_atom)
   use atom_type
   type(atom),intent(inout) :: current_atom

   character(len=128) :: line
   integer :: atomic_number, ios
   character(len=2) :: element
   integer :: radius_single, radius_double, radius_triple
   integer :: unit = 11

100 format(i3,a2,i3,2(i3))
101 format(i3,a2,i3)
102 format(i3,1x,a2,1x,i3,1x,i3,1x,i3)

   open(unit, file="files/CoV_radii", status="old", action="read", iostat=ios)
   if (ios /= 0) then
      print *, "Error opening file!"
      stop
   end if

   do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit  ! EOF
      if (line(1:1) == "#") cycle  ! skip comments

      radius_double = 0
      radius_triple = 0

      read(line, 100, iostat=ios) atomic_number, element, radius_single, radius_double, radius_triple
      if (ios /= 0) then
         read(line, 101, iostat=ios) atomic_number, element, radius_single
         if (ios /= 0) cycle
      end if

      !if the input atom type is the same as the element in the file, update the atom radii
      if (trim(current_atom%element(1:1)) == trim(element)) then
         current_atom%radius1 = radius_single
         current_atom%radius2 = radius_double
         current_atom%radius3 = radius_triple
      end if

      !print 102, atomic_number, element, radius_single, radius_double, radius_triple
   end do

   close(unit)
end subroutine rayon_covalence
