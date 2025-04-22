program ReadCovalentRadii
   implicit none

   ! Define variables
   character(len=256) :: fileName
   character(len=128) :: line
   integer :: atomic_number, ios
   character(len=2) :: element
   integer :: radius_single, radius_double, radius_triple
   integer :: unit = 10

100 format(i3,a2,i3,2(i3))
101 format(i3,a2,i3)
102 format(i3,1x,a2,1x,i3,1x,i3,1x,i3)

   if(iargc() /= 1) then
      print '(a)', "Please provide a file name"
      stop 10
   end if

   call getarg(1, fileName)

   print '(a,a)', "File to read = ", trim(fileName)

   ! Open file
   open(unit, file=fileName, status="old", action="read", iostat=ios)
   if (ios /= 0) then
      print *, "Error opening file!"
      stop
   end if

   do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == "#") cycle  

      radius_double = 0
      radius_triple = 0

      read(line, 100, iostat=ios) atomic_number, element, radius_single, radius_double, radius_triple
      if (ios /= 0) then
         read(line, 101, iostat=ios) atomic_number, element, radius_single
         if (ios /= 0) cycle  ! Skip malformed lines
      end if

      print 102, atomic_number, element, radius_single, radius_double, radius_triple
   end do

   close(unit)

end program ReadCovalentRadii
