module simulation_utils
   use omp_lib
   use atom_type
   use read_mol2
   use molecule_type
   use bond_utils

contains

   function evaluate_fitness(ligand, site) result(score)
      use molecule_type
      use read_mol2
      use bond_utils
      implicit none
      type(molecule), intent(in) :: ligand, site
      real :: score

      score = detect_hydrogen_bonds_between(ligand, ligand%nb_atoms, site, site%nb_atoms)
      !score = 1.0
   end function evaluate_fitness

end module simulation_utils

program genetic_docking
   use simulation_utils
   use io_utils
   implicit none

   integer, parameter :: population_size = 100
   integer, parameter :: max_generations = 500
   real, parameter :: crossover_rate = 0.95
   real, parameter :: mutation_rate = 0.05
   integer, parameter :: log_interval = 5

   integer :: ok

   type(molecule) :: site, ligand
   !type(atom), allocatable :: site(:)
   !type(atom), allocatable :: population(:,:)
   type(molecule), allocatable :: population(:)
   real, allocatable :: fitness(:)
   character(len=128) :: ligand_file, site_file
   integer :: i, gen, site_size, ligand_size
   INTEGER :: maxl(1)
   real :: best_fitness
   !integer :: string_size
   !character(len=256) :: bounds_file
   logical :: with_bond
   real :: hbond_count
   real :: score
   character(len=256), parameter :: output_directory = "results/"
   double precision :: start_time, end_time
   double precision :: ostart,oend

   !legacy argc instrinsic
   integer :: n_args

   n_args = command_argument_count()

   allocate(fitness(population_size), stat=ok)
   if (ok /= 0) then
      print *, "Error during fitness allocation"
      stop 1
   end if
   allocate(population(population_size), stat=ok)
   if (ok /= 0) then
      print *, "Error during population allocation"
      stop 2
   end if
   if (n_args /= 2) then
      print *, "Usage: genetic_docking <ligand_file> <site_file>"
      stop 4
   end if

   ! Read command line arguments
   ! Assuming the first argument is the ligand file and the second is the site file

   call getarg(1, ligand_file)
   call getarg(2, site_file)

   !$OMP PARALLEL
   if(omp_get_thread_num() == 0) then
      print *, "Number of threads = ", omp_get_num_threads()
   end if
   !$OMP END PARALLEL

   call read_molecule_file(ligand_file, ligand, ligand_size)
   call read_molecule_file(site_file, site, site_size)

   call write_bounds_to_file(ligand%atoms, ligand_size, ligand_file)

   print *, "Ligand size: ", ligand_size
   print *, "Site size: ", site_size

   call read_bounds_from_file(site_size, site_file)

   !test bound
   with_bond = is_bonded(site%atoms(1524)%number, site%atoms(1531)%number)
   if (with_bond) then
      print *, "Atoms are bonded"
   else
      print *, "Atoms are not bonded"
   end if

   call initialize_population(ligand, population_size, population)

   hbond_count = detect_hydrogen_bonds_between(ligand, ligand_size, site, site_size)

   !create data log file
   call create_csv("fitness_log.csv")

   !start the genetic algorithm
   call CPU_TIME(start_time) ! Start timing
   ostart = omp_get_wtime()

   do gen = 1, max_generations
      do i = 1, population_size
         score = evaluate_fitness(population(i), site)
         fitness(i) = score
      end do
!
      if (mod(gen, log_interval) == 0) then
         call log_fitness(gen, fitness, population_size)
         !save the best fitness molecule each log_interval generation
         best_fitness = maxval(fitness)
         maxl = maxloc(fitness)
         print *, "Best fitness: ", best_fitness
         print *, "Best index: ", maxl(1)
         call save_best_population(population(maxl(1)), gen, output_directory)
      end if
!
      call evolve_population(population, population_size, fitness, crossover_rate, mutation_rate)
   end do

   call CPU_TIME(end_time) ! End timing
   oend = omp_get_wtime()

   write(*,*) 'Fortran CPU time elapsed', end_time-start_time
   write(*,*) 'OpenMP Walltime elapsed', oend-ostart
end program genetic_docking

subroutine log_fitness(gen, fitness, pop_size)
   use io_utils
   implicit none
   integer, intent(in) :: pop_size
   integer, intent(in) :: gen
   real, intent(in) :: fitness(pop_size)
   integer :: i

   !print *, "Generation: ", gen
   !print *, "Fitness: "
   !do i = 1, size(fitness)
   !   print *, fitness(i)
   !end do
   call add_csv_entry("fitness_log.csv", gen, minval(fitness), maxval(fitness), sum(fitness) / pop_size)
end subroutine log_fitness

subroutine save_best_population(mol, generation, directory)
   use atom_type
   use molecule_type
   use read_mol2
   implicit none

   type(molecule), intent(in) :: mol
   integer, intent(in) :: generation
   character(len=256), intent(in) :: directory
   character(len=256) :: filename

   write(filename, "(A,I0,A)") "best_ligand_", generation, ".mol2"
   call save_ligand_molecule(mol, generation, filename, directory)
end subroutine save_best_population

subroutine initialize_population(base_ligand, pop_size, population)
   use atom_type
   use read_mol2
   implicit none
   integer, intent(in) :: pop_size
   type(molecule), intent(inout) :: population(pop_size)
   type(molecule), intent(in) :: base_ligand !init molecule
   integer :: base_ligand_size
   integer :: i

   base_ligand_size = base_ligand%nb_atoms

   do i = 1, pop_size
      call population(i)%init_mol(base_ligand_size)
      population(i) = base_ligand
      call randomize_conformation(base_ligand, base_ligand_size, population(i))
   end do
end subroutine initialize_population

subroutine randomize_conformation(base, base_size, individual)
   use atom_type
   use molecule_type
   implicit none
   type(molecule), intent(in) :: base
   integer, intent(in) :: base_size
   type(molecule), intent(inout) :: individual
   integer :: i
   real :: dx, dy, dz
   real :: max_displacement
   real :: r

   max_displacement = 0.5  ! max displacement in angstroms

   do i = 1, base_size
      individual%atoms(i)%type = base%atoms(i)%type

      call random_number(r)
      dx = (2.0 * r - 1.0) * max_displacement
      call random_number(r)
      dy = (2.0 * r - 1.0) * max_displacement
      call random_number(r)
      dz = (2.0 * r - 1.0) * max_displacement

      individual%atoms(i)%coordinates(1) = base%atoms(i)%coordinates(1) + dx
      individual%atoms(i)%coordinates(2) = base%atoms(i)%coordinates(2) + dy
      individual%atoms(i)%coordinates(3) = base%atoms(i)%coordinates(3) + dz
   end do
end subroutine randomize_conformation

subroutine evolve_population(pop, pop_size, fitness, crossover_rate, mutation_rate)
   use atom_type
   use molecule_type
   use omp_lib
   implicit none
   !integer, intent(in) :: site_size
   integer, intent(in) :: pop_size
   type(molecule), intent(inout) :: pop(pop_size)
   real, intent(in) :: fitness(pop_size), crossover_rate, mutation_rate
   type(molecule), allocatable :: new_pop(:)
   integer :: n, i, p1, p2, j
   real :: r

   call random_number(r)

   n = size(fitness)
   allocate(new_pop(pop_size))

   do i = 1, n, 2
      call random_number(r)
      call select_parents(fitness, pop_size, p1, p2)
      !print *, "Number of atoms: ", pop(p1)%nb_atoms
      !print *, "Number of atoms: ", pop(p2)%nb_atoms
      if (r < crossover_rate) then
         call crossover(pop(p1), pop(p2), new_pop(i), new_pop(i+1))
      else
         new_pop(i) = pop(p1)
         new_pop(i)%nb_atoms = pop(p1)%nb_atoms
         new_pop(i+1) = pop(p2)
         new_pop(i+1)%nb_atoms = pop(p2)%nb_atoms
      end if
   end do

   !$omp parallel do
   do i = 1, n
      call random_number(r)
      if (r < mutation_rate) then
         call mutate(new_pop(i))
      end if
   end do
   !$omp end parallel do

      !print population
   !do i = 1, n
   !   do j = 1, pop(i)%nb_atoms
   !      print *, "Population: ", i, " ", pop(i)%atoms(j)%element, pop(i)%atoms(j)%coordinates(1), &
   !      pop(i)%atoms(j)%coordinates(2), pop(i)%atoms(j)%coordinates(3)
   !   end do
   !end do

   pop = new_pop
end subroutine evolve_population

subroutine select_parents(fitness, pop_size, p1, p2)
   implicit none
   integer, intent(in) :: pop_size
   real, intent(in) :: fitness(pop_size)
   integer, intent(out) :: p1, p2
   real :: total_fitness, rand_val, cumulative
   integer :: i
   real :: r

   total_fitness = sum(fitness)

   call random_number(r)
   rand_val = r * total_fitness
   cumulative = 0.0
   do i = 1, pop_size
      cumulative = cumulative + fitness(i)
      if (cumulative >= rand_val) then
         p1 = i
         exit
      end if
   end do

   do
      call random_number(r)
      rand_val = r * total_fitness
      cumulative = 0.0
      do i = 1, pop_size
         cumulative = cumulative + fitness(i)
         if (cumulative >= rand_val) then
            if (i /= p1) then
               p2 = i
               exit
            end if
         end if
      end do
      if (p2 /= 0 .and. p2 /= p1) exit
   end do

   !print *, "Parent 1: ", p1
   !print *, "Parent 2: ", p2
end subroutine select_parents

 
subroutine crossover(parent1, parent2, child1, child2)
   use atom_type
   use molecule_type
   implicit none
   !integer, intent(in) :: site_size
   type(molecule), intent(in) :: parent1, parent2
   type(molecule), intent(out) :: child1, child2
   integer :: i, crossover_point
   real :: r

   call random_number(r)

   crossover_point = int(r * parent1%nb_atoms)

   child1%nb_atoms = parent1%nb_atoms
   child2%nb_atoms = parent1%nb_atoms

   allocate(child1%atoms(child1%nb_atoms))
   allocate(child2%atoms(child2%nb_atoms))

   do i = 1, parent1%nb_atoms
      if (i <= crossover_point) then
         child1%atoms(i) = parent1%atoms(i)
         child2%atoms(i) = parent2%atoms(i)
      else
         child1%atoms(i) = parent2%atoms(i)
         child2%atoms(i) = parent1%atoms(i)
      end if
   end do
end subroutine crossover

subroutine mutate(individual)
   use atom_type
   use molecule_type
   implicit none
   !integer, intent(in) :: site_size
   type(molecule), intent(inout) :: individual
   integer :: mutation_point
   real :: r, delta
   real, parameter :: range = 10.0

   call random_number(r)

   mutation_point = r * (individual%nb_atoms - 1.0) + 1.0
   !print *, "Number of atoms: ", individual%nb_atoms
   !print *, "Random number: ", r
   !print *, "Mutation point: ", mutation_point

   ! Randomly change the coordinates of the atom at mutation_point
   call random_number(r)
   delta = (r - 0.5) * range
   !print *, "Delta: ", delta
   
   individual%atoms(mutation_point)%coordinates(1) = individual%atoms(mutation_point)%coordinates(1) - delta
   call random_number(r)
   delta = (r - 0.5) * range

   individual%atoms(mutation_point)%coordinates(2) = individual%atoms(mutation_point)%coordinates(2) - delta
   call random_number(r)
   delta = (r - 0.5) * range

   individual%atoms(mutation_point)%coordinates(3) = individual%atoms(mutation_point)%coordinates(3) - delta
end subroutine mutate
