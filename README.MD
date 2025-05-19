# Projet CHPS1002 - Amerrage/docking en Fortran

Le code est diviser en 2 parties :
* `part1` contenant le code pour la partie du TP *"Les jolis petits rayons"*
* `part2` contenant le code pour la partie du TP *"On danse ensemble ?"*

## Compilation

### Partie 1
```bash
cd part1
make
./complete_mol2_omp ../files/1QSN_NO_BOUND.mol2
```

### Partie 2

```bash
cd part2
make genetic_algorithm
OMP_NUM_THREADS=<nombre_de_threads> ./genetic_algorithm ../files/ligand_NO_BOND.mol2 ../files/site.mol2
```

Les résultats de *docking* de la partie 2 sont disponibles dans le dossier `results`