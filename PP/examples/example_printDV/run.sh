#!/bin/bash

cat >SrVO3.scf.in << EOF
 &control
    calculation='scf'
    prefix = 'SrVO3'
    pseudo_dir = './PP'
    outdir = 'results'
    verbosity = 'high'
 /
 &system
    ibrav = 1
    celldm(1) = 7.226
    nat = 5
    ntyp = 3
    ecutwfc = 100
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
 /
 &electrons
    conv_thr =  1.0d-12
 /
ATOMIC_SPECIES
 Sr  87.62    Sr.upf
 V   50.942   V.upf
 O   15.999   O.upf
ATOMIC_POSITIONS {crystal}
 V   0.000  0.000  0.000
 O   0.500  0.000  0.000
 O   0.000  0.500  0.000
 O   0.000  0.000  0.500
 Sr  0.500  0.500  0.500
K_POINTS {automatic}
 6 6 6   0 0 0
EOF

echo -e "SCF ...\c"
mpirun -np 4 /home/colonna_n/CODES/q-e_for-dyn-Hub/bin/pw.x -in SrVO3.scf.in > SrVO3.scf.out
echo -e " DONE"

cat > SrVO3.bands.in <<EOF 
 &control
    calculation='bands'
    prefix = 'SrVO3'
    pseudo_dir = './PP'
    outdir = 'results'
    verbosity = 'high'
 /
 &system
    ibrav = 1
    celldm(1) = 7.226
    nat = 5
    ntyp = 3
    ecutwfc = 100
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
    nbnd = 40
 /
 &electrons
    conv_thr =  1.0d-12
 /
ATOMIC_SPECIES
 Sr  87.62    Sr.upf
 V   50.942   V.upf
 O   15.999   O.upf
ATOMIC_POSITIONS {crystal}
 V   0.000  0.000  0.000
 O   0.500  0.000  0.000
 O   0.000  0.500  0.000
 O   0.000  0.000  0.500
 Sr  0.500  0.500  0.500
$(sh kmesh.sh 0.1 0.1 0.1)
EOF


echo -e "BANDS ...\c"
mpirun -np 4 /home/colonna_n/CODES/q-e_for-dyn-Hub/bin/pw.x -in SrVO3.bands.in > SrVO3.bands.out
echo -e " DONE"

cat > print_dV.in <<EOF
&pp
   prefix='SrVO3'
   outdir='./results'
/
EOF

echo -e "Print DV ...\c"
/home/colonna_n/CODES/q-e_for-dyn-Hub/bin/print_deltaV.x -in print_dV.in > print_dV.out
echo -e " DONE"

