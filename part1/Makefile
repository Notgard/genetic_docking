CC = gfortran
#CC = ifort

# Compilation flags
FCFLAGS = -c -Jmod -fopenmp
LDFLAGS = -fopenmp

EXEC = complete_mol2_omp molecule atom

complete_mol2_omp: complete_mol2_omp.o molecule_type.o atom_type.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) $(LDFLAGS) -o $@ $+

molecule: molecule_type.o atom_type.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) $(LDFLAGS) -o $@ $+

atom: atom_type.o
	-@echo ""
	-@echo "Linking    $(@)"
	-@echo ""
	-@$(CC) $(LDFLAGS) -o $@ $+

molecule_type.o: atom_type.o
complete_mol2_omp.o: molecule_type.o atom_type.o

%.o: %.f90
	-@echo ""
	-@echo "Generating $@"
	-@$(CC) $(FCFLAGS) $< -o $@

help:
	@echo $(EXEC)
	@echo "(C) G.Roubahie-Fissa"
	@echo "Sous-commandes :"
	@echo "NomFichier (sans son extension)       : compile et fait l'édition de lien du fichier \"NomFichier.f90\" correspondant pour générer un exécutable"
	@echo "make cleanSource                      : supprime les fichiers parasites (*~, *.old,#*,*.bak)"
	@echo "make clean                            : supprime *tous* les fichiers reproductibles ici les .o, les .mod  et aussi les fichiers parasites"
	@echo "make mrproper                         : supprime *tous* les fichiers reproductibles, les fichiers parasites et aussi les exécutables"

###------------------------------
### Cleaning
###------------------------------
clean:
	-@rm -rf *.o *.mod

mrproper: clean cleanSource
	-@rm -rf $(EXEC)

cleanSource:
	-@find . \( -name "*~" -o -name "*.old" -o -name "#*" \) -print -exec rm \{\} \;

.PHONY: $(EXEC) clean mrproper cleanSource
