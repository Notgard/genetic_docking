program LectureFichier
   implicit none

   integer :: i, nbArg, end, ok
   character(len=128) :: fileName, line
   integer :: atomNumber, residueNumber
   character(len=4) :: atomName
   character(len=6) :: atomType
   character(len=6) :: residueName
   real :: x, y, z, charge
   character(len=80) :: mol_format = '(i5, a4, 3f8.4, a6, i5, a6, f8.4)'
   character(len=40) :: print_format = '(a4, 3f8.4)'

   if(iargc() /= 1) then
      print '(a)', "Please provide a file name"
      stop 10
   end if

   call getarg(1, fileName)

   print '(a,a)', "File to read = ", trim(fileName)

   open(unit=10, file=fileName, iostat=ok, status='old')

   if(ok /= 0) then
      print '(a,4x,a)', "Error during opening", fileName
      stop 20
   end if

   do
      read(10, '(A)', iostat=end) line
      if (end /= 0) exit
      if (trim(line) == '@<TRIPOS>ATOM') then
         exit
      end if
   end do

   ! Read ATOM lines
   do
      read(10, *, iostat=end) atomNumber, atomName, x, y, z, atomType, residueNumber, residueName, charge
      if (end /= 0) exit

      print print_format, &
            atomName, x, y, z
   end do

   close(10)

   100 format(i5, a4, 3f8.4, a6, i5, a6, f8.4)
   200 format(a4, 3f8.4)

end program LectureFichier
