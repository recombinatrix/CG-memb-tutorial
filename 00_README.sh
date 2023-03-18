#!/usr/bin/env bash

# Setup replicate of GlyT2 in neuronal membrane

# use go_martinize.py to create GO Martini model of GlyT2

./go_martinize.py -f 01_GlyT2.pdb -o GlyT2_CG.top -x 02_GlyT2_CG.pdb -dssp mkdssp -p backbone -ff elnedyn22 -go -goepsilon 9.414 > go_martinize_GlyT2.out

# attach substrate and ions to 02_GlyT2_CG.pdb

gsed -i '$d' 02_GlyT2_CG.pdb 
gsed -i '$d' 02_GlyT2_CG.pdb 
cat 02_GlyT2_CG.pdb 02_append_sub_ions.pdb > 03_GlyT2_CG_sub_ions.pdb

# edit topology file
# change martini.itp to martini_v2.2.itp
gsed -i 's?martini.itp?martini.ff/martini_v2.2P.itp?' GlyT2_CG.top

# edit topology file to include substrate and ion topology
gsed -i '\?#include "martini.ff/martini_v2.2P.itp"?a #include "martini.ff/martini_v2.0_ions.itp"' GlyT2_CG.top
gsed -i '\?#include "martini.ff/martini_v2.2P.itp"?a #include "martini.ff/glycine-CG.itp"' GlyT2_CG.top

# update molecule list
gsed -i '\?^Protein ?a Glycine  1' GlyT2_CG.top
gsed -i '\?^Glycine ?a NA+   3' GlyT2_CG.top

# generate box
gmx editconf -f 03_GlyT2_CG_sub_ions.pdb -o 04_GlyT2_CG_box.gro -box 20 20 14

# minimise protein in vacuum
gmx grompp -f emin.mdp -c 04_GlyT2_CG_box.gro -p GlyT2_CG.top -o  05_GlyT2_CG_min.tpr -maxwarn 1 2>&1 | tee 05_grompp.out
gmx mdrun -v -deffnm 05_GlyT2_CG_min

# Checking it, it looks ok. The substrate and ions are all present

# center protein in box
printf "1\n0\n" | gmx trjconv -f 05_GlyT2_CG_min.gro  -o 05_GlyT2_CG_min.pdb -center -s  05_GlyT2_CG_min.tpr -pbc mol

# generate membrane
./insane_13July2020.py -f 05_GlyT2_CG_min.pdb -o 06_GlyT2_CG_neuronal.gro -p GlyT2_CG_neuronal.top -d 0 -x 20 -y 20 -z 14 -sol PW -center -l CHOL:0.446 -l PUPE:0.097 -l PAPE:0.061 -l POPC:0.049 -l PUPS:0.034 -l DPPC:0.03 -l PAPS:0.028 -l PAPC:0.026 -l POPE:0.025 -l POPS:0.025 -l PUPI:0.02 -l DPSM:0.015 -l OUPE:0.014 -l OAPE:0.013 -l POPI:0.013 -l PAPI:0.013 -l DOPC:0.012 -l PUPC:0.01 -l OUPS:0.007 -l DPPS:0.005 -l PIPI:0.005 -l DPCE:0.004 -l PFPC:0.003 -l OIPC:0.003 -l OIPE:0.003 -l PNSM:0.003 -l PBSM:0.003 -l PAPA:0.003 -l PAP1:0.003 -l PAP2:0.003 -l PAP3:0.003 -l PADG:0.003 -l OUPC:0.002 -l POSM:0.002 -l POP1:0.002 -l POP2:0.002 -l POP3:0.002 -l IPE:0.002 -l POPA:0.001 -l DBCE:0.001 -l PNCE:0.001 -l PPC:0.001 -l IPC:0.001 -l PPE:0.001 -l PODG:0.001 -u CHOL:0.444 -u POPC:0.087 -u DPSM:0.058 -u DPPC:0.053 -u PUPE:0.05 -u DPGS:0.049 -u PAPC:0.046 -u PAPE:0.031 -u DOPC:0.022 -u PUPC:0.017 -u POPE:0.013 -u PNSM:0.013 -u PBSM:0.011 -u PNGS:0.011 -u DPG1:0.009 -u DPG3:0.009 -u DBGS:0.009 -u OAPE:0.007 -u OUPE:0.007 -u POSM:0.007 -u PFPC:0.006 -u OIPC:0.006 -u POGS:0.006 -u OUPC:0.004 -u DPCE:0.004 -u PADG:0.003 -u DBG1:0.002 -u PNG1:0.002 -u DBG3:0.002 -u PNG3:0.002 -u PPC:0.002 -u OIPE:0.001 -u POG1:0.001 -u POG3:0.001 -u DBCE:0.001 -u PNCE:0.001 -u IPC:0.001 -u PPE:0.001 -u IPE:0.001 -u PODG:0.001

# edit topology file
# add #include lines for itp files, #define rubber bands
gsed -i '\?#include "martini.ff/martini_v2.0_ions.itp"?a #define RUBBER_BANDS\n' GlyT2_CG_neuronal.top
gsed -i '\?#include "martini.ff/martini_v2.0_ions.itp"?a #include "Protein.itp"' GlyT2_CG_neuronal.top
gsed -i '\?#include "martini.ff/martini_v2.0_ions.itp"?a #include "martini.ff/glycine-CG.itp"' GlyT2_CG_neuronal.top

# change sim name
gsed -i 's?Protein in INSANE?GlyT2 in neuronal INSANE?' GlyT2_CG_neuronal.top

# update molecule list
gsed -i '\?^Protein ?a Glycine  1' GlyT2_CG_neuronal.top
gsed -i '\?^Glycine ?a NA+   3' GlyT2_CG_neuronal.top

# minimise protein in membrane
gmx grompp -f emin.mdp  -c 06_GlyT2_CG_neuronal.gro  -p GlyT2_CG_neuronal.top -o 07_GlyT2_neuronal_min.tpr  2>&1 | tee 07_grompp.out
gmx mdrun -v -deffnm 07_GlyT2_neuronal_min


# Neutralize and add 0.15M NaCl
gmx grompp -f emin.mdp -c 07_GlyT2_neuronal_min.gro -p GlyT2_CG_neuronal.top -o 08_GlyT2_neuronal_ions.tpr 2>&1 | tee 08_grompp.out
printf "PW\n" | gmx genion -s 08_GlyT2_neuronal_ions.tpr  -p GlyT2_CG_neuronal.top -conc 0.15 -neutral -o  08_GlyT2_neuronal_ions.gro

# change NA/CL to NA+/CL-
gsed -i 's?NA      NA?ION    NA+?' 08_GlyT2_neuronal_ions.gro
gsed -i 's?CL      CL?ION    CL-?' 08_GlyT2_neuronal_ions.gro
gsed -i 's/^NA /NA+/' GlyT2_CG_neuronal.top
gsed -i 's/^CL /CL-/' GlyT2_CG_neuronal.top

# Minimise again
gmx grompp -f emin.mdp -c 08_GlyT2_neuronal_ions.gro -p GlyT2_CG_neuronal.top -o 09_GlyT2_neuronal_ions_min.tpr 2>&1 | tee 09_grompp.out
gmx mdrun -v -deffnm 09_GlyT2_neuronal_ions_min

# Make index file
printf '!1 & !13 & !63\nname 64 lipids\n13 | 63\nname 65 PW_ION\n\nq\n' | gmx make_ndx -f 09_GlyT2_neuronal_ions_min -o GlyT2_neuronal.ndx

# INDEX GROUPS FOR COUPLING ARE Protein, lipids, and PW_ION
# because my model has ions first due to the GlyT2 sodium ions, gotta rename group 65 from ION_PW to PW_ION for consistency with mdp files.


# Equilibrate
gmx grompp -f 10_eq1000_10fs.mdp -c 09_GlyT2_neuronal_ions_min.gro -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 10_GlyT2_neuronal_eq1000.tpr -r 09_GlyT2_neuronal_ions_min.gro 2>&1 | tee 10_grompp.out
gmx mdrun -v -deffnm 10_GlyT2_neuronal_eq1000

gmx grompp -f 11_eq500_10fs.mdp -c 10_GlyT2_neuronal_eq1000.gro -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 11_GlyT2_neuronal_eq500.tpr -r 10_GlyT2_neuronal_eq1000.gro 2>&1 | tee 11_grompp.out
gmx mdrun -v -deffnm 11_GlyT2_neuronal_eq500

gmx grompp -f 12_eq100_10fs.mdp -c 11_GlyT2_neuronal_eq500.gro -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 12_GlyT2_neuronal_eq100.tpr -r 11_GlyT2_neuronal_eq500.gro 2>&1 | tee 12_grompp.out
gmx mdrun -v -deffnm 12_GlyT2_neuronal_eq100

gmx grompp -f 13_eq50_10fs.mdp -c 12_GlyT2_neuronal_eq100.gro -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 13_GlyT2_neuronal_eq50.tpr -r 12_GlyT2_neuronal_eq100.gro 2>&1 | tee 13_grompp.out
gmx mdrun -v -deffnm 13_GlyT2_neuronal_eq50

gmx grompp -f 14_eq10_10fs.mdp -c 13_GlyT2_neuronal_eq50.gro -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 14_GlyT2_neuronal_eq10.tpr -r 13_GlyT2_neuronal_eq50.gro 2>&1 | tee 14_grompp.out
gmx mdrun -v -deffnm 14_GlyT2_neuronal_eq10

gmx grompp -f 15_eq0_10fs.mdp -c 14_GlyT2_neuronal_eq10.gro -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 15_GlyT2_neuronal_eq0_10fs.tpr -r 14_GlyT2_neuronal_eq10.gro 2>&1 | tee 15_grompp.out
gmx mdrun -v -deffnm 15_GlyT2_neuronal_eq0_10fs

gmx grompp -f 16_eq0_15fs.mdp -c 15_GlyT2_neuronal_eq0_10fs -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 16_GlyT2_neuronal_eq0_15fs.tpr -r 15_GlyT2_neuronal_eq0_10fs 2>&1 | tee 16_grompp.out
gmx mdrun -v -deffnm 16_GlyT2_neuronal_eq0_15fs

gmx grompp -f 17_eq0_20fs.mdp -c 16_GlyT2_neuronal_eq0_15fs -n GlyT2_neuronal.ndx -p GlyT2_CG_neuronal.top -o 17_GlyT2_neuronal_eq0_20fs.tpr -r 16_GlyT2_neuronal_eq0_15fs -maxwarn 1  2>&1 | tee 17_grompp.out
gmx mdrun -v -deffnm 17_GlyT2_neuronal_eq0_20fs

# Step 17  grompp (20 fs timestep unrestrained EQ) has been failing because of the following warning 
#
#     WARNING 1 [file GlyT2_CG_neuronal.top, line 114]:
#      The bond in molecule-type Protein between atoms 1 BB and 4 BB has an
#      estimated oscillational period of 9.7e-02 ps, which is less than 5 times
#      the time step of 2.0e-02 ps.
#      Maybe you forgot to change the constraints mdp option.
#
# Consequently I'm running step 17 with a -maxwarn 1 flag, but that means you need to check the grompp output to see what warnings were generated, and check the EQ simulation was reasonable

echo "Final unrestrained EQ was run with the -maxwarn 1 flag, make sure you check step 17 is sensible"
