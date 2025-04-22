module bond_utils
   implicit none
   logical, allocatable :: bond_matrix(:,:)
contains

   function is_bonded(atom1, atom2) result(bond_exists)
      integer, intent(in) :: atom1, atom2
      logical :: bond_exists
      if (.not. allocated(bond_matrix)) then
         bond_exists = .false.
      else
         bond_exists = bond_matrix(atom1, atom2)
      end if
   end function is_bonded
 
end module bond_utils

module read_mol2
   use atom_type
   use molecule_type
   use omp_lib 
contains
   subroutine read_molecule_file(fileName, mol, atom_nums)
      implicit none

      character(len=128), intent(in) :: fileName
      type(molecule), intent(inout) :: mol 
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
      integer, parameter :: NB_MAX_ATOMS=3000

      ! Initialize atoms and counter
      !if (allocated(atoms)) deallocate(atoms)
      !allocate(atoms(0))
      atom_nums = 0

      ! Initialize molecule
      call mol%init_mol(NB_MAX_ATOMS)

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
         call mol%add_atom(a)
         !atoms = [atoms, a]
         atom_nums = atom_nums + 1
      end do

      close(unit)
   end subroutine read_molecule_file

   subroutine rayon_covalence(current_atom)
      use atom_type
      type(atom),intent(inout) :: current_atom

      character(len=128) :: line
      integer :: atomic_number, ios
      character(len=2) :: element
      integer :: radius_single, radius_double, radius_triple
      integer :: unit = 11

100   format(i3,a2,i3,2(i3))
101   format(i3,a2,i3)
!102   format(i3,1x,a2,1x,i3,1x,i3,1x,i3)

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

103   format(1x, i5, 1x, i5, 1x, i5, 1x, i1)

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

   subroutine read_bounds_from_file(atom_nums, fileName)
      use atom_type
      use bond_utils 
      implicit none

      integer, intent(in) :: atom_nums
      character(len=128), intent(in) :: fileName

      integer :: unit = 14
      integer :: ios
      character(len=256) :: line
      integer :: bond_id, atom1, atom2
      character(len=3) :: bond_type
      integer :: num_bonds = 0
      logical :: is_bond_section = .false.

      print '(a,a)', "Reading bonds from file: ", trim(fileName)
      open(unit=unit, file=fileName, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         print *, "Error opening file for bond reading: ", trim(fileName)
         return
      end if

      ! Allocate and initialize bond matrix
      if (.not. allocated(bond_matrix)) then
         allocate(bond_matrix(atom_nums, atom_nums))
         bond_matrix = .false.
      end if

      ! Find bond section
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (trim(line) == '@<TRIPOS>BOND') then
            is_bond_section = .true.
            exit
         end if
      end do

      ! Read bonds if section found
      if (is_bond_section) then
         do
            read(unit, *, iostat=ios) bond_id, atom1, atom2, bond_type
            if (ios /= 0) exit

            num_bonds = num_bonds + 1
            if (atom1 <= atom_nums .and. atom2 <= atom_nums) then
               bond_matrix(atom1, atom2) = .true.
               bond_matrix(atom2, atom1) = .true.
            else
               print *, "Warning: Bond references invalid atom indices:", atom1, atom2
            end if
         end do
      else
         print *, "No bond section found in file:", trim(fileName)
      end if

      close(unit)
      print *, "Read", num_bonds, "bonds from file"
   end subroutine read_bounds_from_file

   !Pour déterminer un pont hydrogène, il faut vérifier qu'un atome d'hydrogène d'un des deux systèmes est à une distance comprise dans l'intervalle [2,2;4] angstroms avec un atome électronégatif et que l'angle formé par cet liaison est dans l'intervalle [90,150] degrés. Le mesure de l'angle se fait selon 3 éléments : l'atome d'hydrogène, l'atome électronégatif et un troisième atome lié à l'atome électronégatif.
   real function detect_hydrogen_bonds_between(ligand, n_lig, site, n_site) result(hbond_count)
      use atom_type
      use molecule_type
      use bond_utils
      use omp_lib
      implicit none

      !type(atom), dimension(:), intent(in) :: ligand, site
      type(molecule), intent(in) :: ligand, site
      integer, intent(in) :: n_lig, n_site

      integer :: i, j, k
      real :: distance, angle_deg
      logical :: electronegs_lig(n_lig), hydrogens_lig(n_lig)
      logical :: electronegs_site(n_site), hydrogens_site(n_site)
      real :: x1, y1, z1, x2, y2, z2, x3, y3, z3
      real :: vec1(3), vec2(3), vec_len1, vec_len2, dot_product
      logical :: with_bond

      hbond_count = 0

      ! Identify electronegative and hydrogen atoms
      electronegs_lig = .false.
      electronegs_site = .false.
      hydrogens_lig = .false.
      hydrogens_site = .false.

      do i = 1, n_lig
         if (trim(ligand%atoms(i)%element(1:1)) == "H") hydrogens_lig(i) = .true.
         if (trim(ligand%atoms(i)%element(1:1)) == "N" .or. trim(ligand%atoms(i)%element(1:1)) == "O" .or. &
            trim(ligand%atoms(i)%element(1:1)) == "F" .or. trim(ligand%atoms(i)%element(1:2)) == "Cl") electronegs_lig(i) = .true.
      end do

      do i = 1, n_site
         if (trim(site%atoms(i)%element(1:1)) == "H") hydrogens_site(i) = .true.
         if (trim(site%atoms(i)%element(1:1)) == "N" .or. trim(site%atoms(i)%element(1:1)) == "O" .or. &
            trim(site%atoms(i)%element(1:1)) == "F" .or. trim(site%atoms(i)%element(1:2)) == "Cl") electronegs_site(i) = .true.
      end do

      ! Main loop: ligand H to site EN, and site H to ligand EN
      !$omp parallel default(shared) private(i,j,k,distance,x1,y1,z1,x2,y2,z2,x3,y3,z3), &
      !$omp private(vec1,vec2,vec_len1,vec_len2,dot_product,angle_deg) reduction(+:hbond_count)
      ! Case 1: H from ligand, EN from site
      do i = 1, n_lig
         if (hydrogens_lig(i)) then
            x1 = ligand%atoms(i)%coordinates(1)
            y1 = ligand%atoms(i)%coordinates(2)
            z1 = ligand%atoms(i)%coordinates(3)

            !print *, "Ligand H: ", ligand%atoms(i)%number, x1, y1, z1

            do j = 1, n_site
               if (electronegs_site(j)) then
                  x2 = site%atoms(j)%coordinates(1)
                  y2 = site%atoms(j)%coordinates(2)
                  z2 = site%atoms(j)%coordinates(3)
                  distance = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

                  !print *, "Distance: ", distance
                  if (distance >= 2.2 .and. distance <= 4.0) then
                     do k = 1, n_site
                        if (k /= j) then
                           x3 = site%atoms(k)%coordinates(1)
                           y3 = site%atoms(k)%coordinates(2)
                           z3 = site%atoms(k)%coordinates(3)
                           !detect bound using function
                           with_bond = is_bonded(site%atoms(j)%number, site%atoms(k)%number)
                           !print *, "Finding Bonded: ", site%atoms(j)%number, site%atoms(k)%number
                           if (with_bond) then
                              !print *, "Bonded: ", site%atoms(j)%number, site%atoms(k)%number
                              ! Calculate angle
                              vec1 = (/x1-x2, y1-y2, z1-z2/)
                              vec2 = (/x3-x2, y3-y2, z3-z2/)
                              vec_len1 = normalization(vec1)
                              vec_len2 = normalization(vec2)
                              if (vec_len1 > 0 .and. vec_len2 > 0) then
                                 dot_product = dot_product3(vec1, vec2)
                                 angle_deg = acos(dot_product / (vec_len1*vec_len2)) * 180.0 / 3.14159265
                                 if (angle_deg >= 90.0 .and. angle_deg <= 150.0) then
                                    !$omp critical
                                       !print *, "H-bond found: ", ligand%atoms(i)%number, site%atoms(j)%number, &
                                       !site%atoms(k)%number, distance, angle_deg
                                    !$omp end critical
                                    hbond_count = hbond_count + 1
                                    exit
                                 end if
                              end if
                           end if
                        end if
                     end do
                  end if
               end if
            end do
         end if
      end do

      ! (Optional) Case 2: H from site, EN from ligand (mirror the loop above)

      !$omp end parallel

      !print *, "Total hydrogen bonds found: ", hbond_count
   end function detect_hydrogen_bonds_between

   ! --- Helpers (à mettre dans un module si besoin)
   pure function normalization(v) result(norm)
      real, dimension(3), intent(in) :: v
      real :: norm
      norm = sqrt(sum(v**2))
   end function normalization

   pure function dot_product3(v1, v2) result(dot)
      real, dimension(3), intent(in) :: v1, v2
      real :: dot
      dot = sum(v1 * v2)
   end function dot_product3

end module read_mol2
