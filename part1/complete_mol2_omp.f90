module molecule_io
   use atom_type
   use omp_lib 
   implicit none

   contains
   
   subroutine read_molecule_file(fileName, atoms, atom_nums)
      implicit none
      
      character(len=128), intent(in) :: fileName
      type(atom), dimension(:), allocatable, intent(out) :: atoms
      integer, intent(out) :: atom_nums
      
      integer :: unit = 10
      integer :: end, ok
      character(len=128) :: line
      integer :: atomNumber, subStructureId
      character(len=4) :: atomName
      character(len=6) :: atomType
      character(len=6) :: subStructureName
      real :: x, y, z, charge
      real, dimension(3) :: coordinates
      type(atom) :: a
      
      ! Initialize atoms and counter
      if (allocated(atoms)) deallocate(atoms)
      allocate(atoms(0))
      atom_nums = 0
      
      print '(a,a)', "File to read = ", trim(fileName)
      open(unit=unit, file=fileName, iostat=ok, status='old')
      if(ok /= 0) then
          print '(a,4x,a)', "Error during opening", fileName
          stop 20
      end if
      
      ! Find atom section
      do
          read(unit, '(a)', iostat=end) line
          if (end /= 0) exit
          if (trim(line) == '@<TRIPOS>ATOM') then
              exit
          end if
      end do
      
      ! Read atoms
      do
          read(unit, *, iostat=end) atomNumber, atomName, x, y, z, atomType, subStructureId, subStructureName, charge
          if (end /= 0) exit
          
          coordinates = (/x, y, z/)
          call a%init_atom(atomNumber, atomName, coordinates, atomType)
          call rayon_covalence(a)
          atoms = [atoms, a]
          atom_nums = atom_nums + 1
      end do
      
      close(unit)
   end subroutine read_molecule_file
   
   ! You can also move write_bounds_to_file here if desired
end module molecule_io

program complete_mol2
   use omp_lib
   use atom_type
   use molecule_io
   implicit none
   integer :: i, j, nbArg, end, ok
   character(len=128) :: fileName, line
   integer :: atomNumber, subStructureId
   character(len=4) :: atomName
   character(len=6) :: atomType
   character(len=6) :: subStructureName
   real :: x, y, z, xb, yb, zb, charge
   real, dimension(3) :: coordinates
   character(len=80) :: mol_format = '(i5, a4, 3f8.4, a6, i5, a6, f8.4)'
   double precision :: start_time, end_time
   double precision :: ostart,oend
   integer :: atom_nums = 0
   integer :: unit = 10
   type(atom) :: a
   type(atom), dimension(:), allocatable :: atoms
   allocate(atoms(0))

   !$OMP PARALLEL
   if(omp_get_thread_num() == 0) then
      print *, "Number of threads = ", omp_get_num_threads()
   end if
   !$OMP END PARALLEL

   if(iargc() /= 1) then
      print '(a)', "Please provide a file name"
      stop unit
   end if

   call getarg(1, fileName)

   call CPU_TIME(start_time) ! Start timing
   ostart = omp_get_wtime()

   call read_molecule_file(fileName, atoms, atom_nums)

   call write_bounds_to_file(atoms, atom_nums, fileName)

   call CPU_TIME(end_time) ! End timing
   oend = omp_get_wtime()

   write(*,*) 'Fortran CPU time elapsed', end_time-start_time
   write(*,*) 'OpenMP Walltime elapsed', oend-ostart

   deallocate(atoms)
end program complete_mol2

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

   open(unit, file="../files/CoV_radii", status="old", action="read", iostat=ios)
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

      !if the input atom type is the same as the element in the file, update the atom radii in angstroms
      if (trim(current_atom%element(1:1)) == trim(element)) then
         if(radius_single == 0) then
            current_atom%radius1 = radius_single
         else
            current_atom%radius1 = radius_single
            current_atom%radius1 = current_atom%radius1 / 100
         end if
         if(radius_double == 0) then
            current_atom%radius2 = radius_double
         else
            current_atom%radius2 = radius_double
            current_atom%radius2 = current_atom%radius2 / 100
         end if
         if(radius_triple == 0) then
            current_atom%radius3 = radius_triple
         else
            current_atom%radius3 = radius_triple
            current_atom%radius3 = current_atom%radius3 / 100
         end if
      end if
      !print 102, atomic_number, element, radius_single, radius_double, radius_triple
   end do

   close(unit)
end subroutine rayon_covalence

subroutine write_bounds_to_file(atoms, atom_nums, fileName)
   use atom_type
   use omp_lib
   implicit none

   type(atom), dimension(atom_nums), intent(in)  :: atoms
   integer, intent(in) :: atom_nums
   character(len=128), intent(in) :: fileName

   integer :: i, j, ios
   real :: radius1_sum, radius2_sum, radius3_sum
   real :: delta = 0.1
   real :: distance, lb, ub
   integer :: old_file_unit = 11
   integer :: new_file_unit = 12
   real :: x, y, z, xb, yb, zb
   integer :: bounds_count = 1
   character(len=256) :: newFileName, line
   integer :: ssize

103 format(1x, i5, 1x, i5, 1x, i5, 1x, i1)

   ssize = len_trim(fileName)
   newFileName = trim(fileName(1:ssize-5)) // "_WITH_BONDS.mol2"
   print *, "Writing to file: ", trim(newFileName)

   open(unit=new_file_unit, file=trim(newFileName), status="replace", action="write", iostat=ios)
   if (ios /= 0) then
      print *, "Error opening file: ", trim(newFileName)
      return
   end if

   !old file to append into the new file
   open(unit=old_file_unit, file=trim(fileName), status="old", action="read", iostat=ios)
   if (ios /= 0) then
      print *, "Error opening original file: ", trim(fileName)
      return
   end if

   !copy contents from original into new file
   do
      read(old_file_unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      write(new_file_unit, '(A)') trim(line)
   end do

   !close original file
   close(old_file_unit)

   !add @<TRIPOS>BOND to the file
   write(new_file_unit, '(A)') "@<TRIPOS>BOND"

   !$omp parallel private(j, x, y, z, xb, yb, zb, distance, radius1_sum, radius2_sum, radius3_sum, lb, ub)
   !$omp do schedule(dynamic)
   do i = 1, atom_nums
      do j = i+1, atom_nums
         x = atoms(i)%coordinates(1)
         y = atoms(i)%coordinates(2)
         z = atoms(i)%coordinates(3)
         xb = atoms(j)%coordinates(1)
         yb = atoms(j)%coordinates(2)
         zb = atoms(j)%coordinates(3)

         distance = sqrt((x - xb)**2 + (y - yb)**2 + (z - zb)**2)
         radius1_sum = atoms(i)%radius1 + atoms(j)%radius1
         radius2_sum = atoms(i)%radius2 + atoms(j)%radius2
         radius3_sum = atoms(i)%radius3 + atoms(j)%radius3

         lb = radius1_sum * (1 - delta)
         ub = radius1_sum * (1 + delta)
         if(distance >= lb .and. distance <= ub) then
            !$omp critical
            write(new_file_unit, 103) bounds_count, atoms(i)%number, atoms(j)%number, 1
            bounds_count = bounds_count + 1
            !$omp end critical
         end if

         lb = radius2_sum * (1 - delta)
         ub = radius2_sum * (1 + delta)
         if(distance >= lb .and. distance <= ub) then
            !$omp critical
            write(new_file_unit, 103) bounds_count, atoms(i)%number, atoms(j)%number, 2
            bounds_count = bounds_count + 1
            !$omp end critical
         end if

         lb = radius3_sum * (1 - delta)
         ub = radius3_sum * (1 + delta)
         if(distance >= lb .and. distance <= ub) then
            !$omp critical
            write(new_file_unit, 103) bounds_count, atoms(i)%number, atoms(j)%number, 3
            bounds_count = bounds_count + 1
            !$omp end critical
         end if
      end do
   end do
   !$omp end do
   !$omp end parallel
   close(new_file_unit)
end subroutine write_bounds_to_file