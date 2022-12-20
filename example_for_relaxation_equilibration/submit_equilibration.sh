#!/bin/bash
#$ -pe mpi 8
#$ -cwd
#$ -N C2EO4_2022-12-19_09_12_06
#$ -j y
#$ -o run.log
#$ -m e
#$ -M test
# requesting 300 hrs wall clock time
#$ -l h_rt=300:00:00
module load openmpi/4.1.1
mpirun -np 8 ~/lammps-29-oct-2020/src/lmp_mpi -in equilibration_run.in > output.txt