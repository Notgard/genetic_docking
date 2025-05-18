# Compiler settings
FC = gfortran
#FC = ifort
#FC = flang

# Compilation flags
FCFLAGS_DEBUG = -c -pg -g -fbacktrace -fbounds-check -fcheck=all -Wall -O0 -Jmod -Imod -fopenmp
FCFLAGS = -O3 -c -Jmod -Imod -fopenmp
LDFLAGS = -fopenmp -pg -g

FCFLAGS = $(FCFLAGS_DEBUG)

# Executables
EXEC = complete_mol2_omp_mod molecule atom genetic_algorithm

# Default target
all: $(EXEC)

# Binaries
complete_mol2_omp_mod: complete_mol2_omp_mod.o molecule_type.o atom_type.o
	@echo ""
	@echo "Linking    $(@)"
	@echo ""
	$(FC) $(LDFLAGS) -o $@ $^

molecule: molecule_type.o atom_type.o
	@echo ""
	@echo "Linking    $(@)"
	@echo ""
	$(FC) $(LDFLAGS) -o $@ $^

io_utils: io_utils.o
	@echo ""
	@echo "Linking    $(@)"
	@echo ""
	$(FC) $(LDFLAGS) -o $@ $^

atom: atom_type.o
	@echo ""
	@echo "Linking    $(@)"
	@echo ""
	$(FC) $(LDFLAGS) -o $@ $^

genetic_algorithm: genetic_algorithm.o complete_mol2_omp_mod.o io_utils.o molecule_type.o atom_type.o
	@echo ""
	@echo "Linking    $(@)"
	@echo ""
	$(FC) $(LDFLAGS) -o $@ $^

# Object files with dependencies
molecule_type.o: molecule_type.f90 atom_type.o
	@echo ""
	@echo "Generating $@"
	$(FC) $(FCFLAGS) $< -o $@

atom_type.o: atom_type.f90
	@echo ""
	@echo "Generating $@"
	$(FC) $(FCFLAGS) $< -o $@

io_utils.o: io_utils.f90
	@echo ""
	@echo "Generating $@"
	$(FC) $(FCFLAGS) $< -o $@

complete_mol2_omp_mod.o: complete_mol2_omp_mod.f90 molecule_type.o atom_type.o
	@echo ""
	@echo "Generating $@"
	$(FC) $(FCFLAGS) $< -o $@

genetic_algorithm.o: genetic_algorithm.f90 complete_mol2_omp_mod.o molecule_type.o atom_type.o
	@echo ""
	@echo "Generating $@"
	$(FC) $(FCFLAGS) $< -o $@

# Help
help:
	@echo $(EXEC)
	@echo "(C) G.Roubahie-Fissa"
	@echo "Sous-commandes :"
	@echo "make                                   : compile tous les exécutables"
	@echo "make <executable>                      : compile seulement l'exécutable spécifié"
	@echo "make cleanSource                       : supprime les fichiers parasites (*~, *.old,#*,*.bak)"
	@echo "make clean                             : supprime *tous* les fichiers reproductibles ici les .o, les .mod et aussi les fichiers parasites"
	@echo "make mrproper                          : supprime *tous* les fichiers reproductibles, les fichiers parasites et aussi les exécutables"

# Cleaning
clean:
	@rm -rf *.o *.mod

mrproper: clean cleanSource
	@rm -rf $(EXEC)

cleanSource:
	@find . \( -name "*~" -o -name "*.old" -o -name "#*" -o -name "*.bak" \) -print -exec rm {} \;

.PHONY: all help clean mrproper cleanSource $(EXEC)
