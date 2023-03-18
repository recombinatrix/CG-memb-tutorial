# CG_memb_tutorial


This tutorial written by Ada Quinn, a researcher in computational chemisty at The University of Queensland.  This readme is version 1.1.0, dated 18 March 2023.  It is distributed under a GNU GPL-v3 (2007) licence.  It is adapted from methods used in a 2021 paper by Katie Wilson, Lily Wang, YC Lin and Megan O'Mara: [Investigating the lipid fingerprint of SLC6 neurotransmitter transporters: a comparison of dDAT, hDAT, hSERT, and GlyT2](https:/doi.org/10.1016/j.bbadva.2021.100010)

I wrote this to provide to students in my own research group, where I would be on hand to support them if something went wrong. If you are not my student please know that while I hope you find it helpful, you use it at your own risk.

## Introduction

This is a tutorial for building a simulation system to study the human glycine transporter GlyT2 (uniprot id [SLC6A5](https://www.uniprot.org/uniprotkb/Q9Y345/entry)) in a complex neuronal membrane.  It was written to work with Gromacs version 2021.4, using the [Martini 2.2P](http://cgmartini.nl/) coarse grained forcefield.

In this tutorial, you will convert a homology model of GlyT2 to a coarse grained representation using [*go_martinize.py*](https://doi.org/10.1021/acs.jctc.6b00986), and then embed this into a binary POPC/cholesterol membrane using a modified version of the [insane.py script](https://doi.org/10.1021/acs.jctc.5b00209).  You will then introduce ions ar physiologically relevant concentration, and perform a series of equilibration steps.  

This entire tutorial and associated files is availible for download from github at [https://github.com/recombinatrix/CG_memb_tutorial](https://github.com/recombinatrix/CG_memb_tutorial).  If you have git installed on your computer, you can download the tutorial and all associated files by running

~~~s
git clone git@github.com:recombinatrix/CG_memb_tutorial.git
~~~

### Requirements

This tutorial assumes you are working with linux and have basic linux knowledge.

You will need working installations of Gromacs 2021.4, Visualising Molecular Dynamics (VMD), and python 2.

You are assumed to know the basics of using VMD.  If you need to learn how to use VMD, I recommend working through [these tutorials](https://www.ks.uiuc.edu/Training/Tutorials/vmd-index.html) from the [Theoretical and Computational Biophysics Group at the University of Illinois Urbana-Champaign](https://www.ks.uiuc.edu/)

I will assume you have already worked through my [united atom membrane system tutorial](https://github.com/recombinatrix/UA-memb-tutorial).

### Useful resources

The [Gromacs manual](https://manual.gromacs.org/2021.4/index.html) includes extensive documentation for every gromacs command used in this tutorial.  Whenever you use a new gromacs command, you should have a look at the manual to see what it does, and see if you can figure out why I wrote it in that particular way.

I strongly encourage you to learn to work with an integrated development environment, such as [Visual Studio Code](https://code.visualstudio.com/).  IDEs allow you to explore your file system, edit multiple documents at once, and run multiple terminal sessions all in a single window.  There are helpful extensions that make allow vscode to understand the shape of gromacs file types, which makes them significantly easier to read and understand.  This is especially valuable while you are learning how to work with gromacs and the associated file formats for the first time.

I suggest that you install VSCode and a [gromacs helper extension](https://marketplace.visualstudio.com/items?itemName=SupernovaZJX.gmx-helper) before you start this tutorial.

### Recordkeeping

Whenever you use a computer to perform scientific process, such as building or analysing a molecular dynamics simulation, you ***must*** make a readme file.  A readme is a written record of everything you did, including every command you ran and every file you changed, along with a minimal explanation as to what you were trying to accomplish.  

Your readme file should be sufficiently detailed that another scientists (who does not know you and does not know what you intended to do) could sit down at your computer and replicate your work with only your readme file.  Make sure your readme is good, because sometimes that scientist who does not know you, and does not know what you intended to do, will be you in the future, trying to figure out what you did six months ago!  

You need to make your own readme file while you work through this tutortial.  **Do not try to write your readme file after you have finished.  That way lies misery.**

## The Tutorial
### Getting started

Any simulation needs three things:

* A set of *instructions*, which tells gromacs how to perform the simulation.  These are usually contained in `.mdp` and `.sh` files, and inside the gromacs code.
* A set of *coordinates*, which describes the position of every atom in the system.  For the kind of work I do, these are usually contained in a `.pdb` or `.gro` file.
* A *forcefield* and associated *parameters*, which describes the forces between atoms, and how different atoms are connected.   These are usually contained in a `.top` file, a number of `.itp` files, and a large number of files contained in a `.ff` folder.

For your tutorial, some of these things are provided.

The *instructions* for running each of the molecular dynamics steps have been provided through a series of `.mdp` and `.sh` files.

The *coordinates* for an atomistic model of GlyT2 have been provided for you.  Over the course of this tutorial you will need to convert them to a coarse grained representation, and use those to generate a system with GlyT2 embedded in a complex membrane.

`01_GlyT2.pdb` contains coordinates of a homology model of GlyT2 based on the homologue dDAT, [published here](https://doi.org/10.1371/journal.pone.0157583).  `02_append_sub_ions.pdb` contains coarse grained coordinates for zwitterionic glycine substrate and bound sodium ions positioned to sit at their respective bindign sites in GlyT2.

The *forcefield* you will use is the [Martini 2.2P](http://cgmartini.nl/) forcefield, which has been provided in the `martini.ff/` folder, along with some useful `.itp` files for various coarse grained lipid species. We will also need to generate *parameters* for additional molecules we introduce into our system, like our coarse-grained representation of GlyT2.

### Use `go_martinize.py` to create GO Martini model of GlyT2

[`go_martinize.py` can be downloaded from the developers](http://info.ifpan.edu.pl/~panos/panos/GoMartini.html), who provide an excellent tutorial on the program.  For more details you should read the [*go_martinize* paper](https://doi.org/10.1021/acs.jctc.6b00986).  You may need to use `chmod` to ensure `go_martinize.py` is executable.  

When building your own systems in future work, you may wish to consider the updated script [`go_martinize_MLO.py`](https://github.com/OkazakiLab/Go-MARTINI), which uses a new cutoff scheme to define native contacts of proteins [Mahmood, Poma and Okazaki, Front. Mol. Biosci. (2021)](https://doi.org/10.3389/fmolb.2021.619381).

`go_martinize.py` uses `dssp`, which is redistributed here under a BSD-2-Clause License.  If you are running this on an apple silicon mac, you may wish to consider [mkdssp](https://github.com/AlexTaguchi/dssp-notes) instead.  As this involves compiling, do so at your own risk.

Use `go_martinize.py` to generate coarse-grained GlyT2.

~~~s
python2 go_martinize.py -f 01_GlyT2.pdb -o GlyT2_CG.top -x      \ 
        02_GlyT2_CG.pdb -dssp ./dssp -p backbone -ff elnedyn22  \
        -go -goepsilon 9.414 > go_martinize_GlyT2.out
~~~

`go_martinize.py` has created several files.  

 * `02_GlyT2_CG.pdb` contained coordinates for coarse grained GlyT2
 * `GlyT2_CG.top` contains the topology of your system
 * `Protein.itp` contains parameters for GlyT2.  If you scroll down to the bottom, you will see a set of position restraints defined by a variable, similar to those you made in my united atom tutorial.

Next we want to place the substrate and ions inside coarse-grained GlyT2 `02_GlyT2_CG.pdb`

~~~s
cat 02_GlyT2_CG.pdb 02_append_sub_ions.pdb > 03_GlyT2_CG_sub_ions.pdb
~~~

Open `03_GlyT2_CG_sub_ions.pdb` and scroll to the end of the document.

~~~s
ATOM   1209   BB PRO   756      67.760  92.360  43.410  1.00  0.00
ATOM   1210  SC1 PRO   756      66.153  92.690  42.520  1.00  0.00
ATOM   1211   BB ASP   757      71.520  93.060  43.590  1.00  0.00
ATOM   1212  SC1 ASP   757      73.531  92.461  44.984  1.00  0.00
TER
ENDMDL
ATOM      1  BB  GLY   759      77.712  84.279  75.577  1.00  0.00           B 
ATOM      2  NA+ ION   760      74.820  83.510  76.000  1.00  0.00            
ATOM      3  NA+ ION   761      79.240  81.380  74.090  1.00  0.00            
ATOM      4  NA+ ION   762      68.150  94.120  60.500  1.00  0.00            
TER
ENDMDL
~~~

Remove the `TER` and `ENDMDL` lines from the middle of this section.

We want to un an energy minimization, but first we need to edit the topology file.

We need to change the forcefield to the Martini 2.2P forcefield. Change the line

~~~c
    #include "martini.itp"
~~~

to

~~~c
    #include "martini.ff/martini_v2.2P.itp"
~~~

You also need to add the parameters for the substrate and ions.  Write inlude statements for the files `martini.ff/martini_v2.0_ions.itp` and `martini.ff/glycine-CG.itp`, and add a glycine molecule and three sodium ions to the `[ molecules ]` list.

Next, use `editconf` to generate a boxm and minimise our goarse grained protein in vaccum

~~~s
gmx editconf -f 03_GlyT2_CG_sub_ions.pdb -o 04_GlyT2_CG_box.gro \
               -box 20 20 14

gmx grompp -f emin.mdp -c 04_GlyT2_CG_box.gro -p GlyT2_CG.top  \
               -o  05_GlyT2_CG_min.tpr 2>&1 | tee 05_grompp.out

gmx mdrun -v -deffnm 05_GlyT2_CG_min
~~~

Look at it in VMD.  Does it look ok?  If it does, use trjconv to center GlyT2 in the box.

~~~s
gmx trjconv -f 05_GlyT2_CG_min.gro  -o 05_GlyT2_CG_min.pdb -center -s  05_GlyT2_CG_min.tpr -pbc mol
~~~

Chose to center the protein, and export the entire system.

### Build a system with the protein in a membrane

We're going to generate membrane with [insane.py script](https://doi.org/10.1021/acs.jctc.5b00209).  We're using a version of insane.py that has been edited to include topologies for some additional lipids.  Insane will build us a membrane, based on a list of lipids we specify in the upper and lower leaflets.  The position of the various lipids will be randomized.  To take advantage of this randomness, when you build a real system you need to repeat this entir process for every replicate.

The full command we're going to use to build this system is enourmous, so we're going to look at a simplified version first.  Run this command to display the insane.py help

~~~s
python2 insane_13July2020.py -h
~~~

Use this to try to understand all the different pieces of the next command.  Don't go ahead until you think you know what all these different flags are doing!

First we're going to make an example system, with GlyT2 in a simple binary membrane.

~~~s
python2 ./insane_13July2020.py -f 05_GlyT2_CG_min.pdb                       \ 
      -o 06_GlyT2_CG_example.gro -p GlyT2_CG_example.top                    \ 
      -d 0 -x 20 -y 20 -z 14 -sol PW -center -l CHOL:0.2 -l POPC:0.8        \
      -u CHOL:0.1 -u POPC:0.8
~~~

Have a look at `06_GlyT2_CG_example.gro` in VMD.  You should be able to see your protein sitting in the middle of a bilayer of POPC and cholesterol, surrounded by coarse grain polarizable water (`resname PW`).  Try editing the different flags in your `insane.py` command and seeing how it changes the output files.

Once you are satisfied you understand what all of the different flags are doing, you're ready to build your real system, in a complex neuronal membrane.  This command is very similar to the example system, but it has so many more lipid species.

~~~s
python2 ./insane_13July2020.py -f 05_GlyT2_CG_min.pdb                       \ 
      -o 06_GlyT2_CG_neuronal.gro -p GlyT2_CG_neuronal.top        \ 
      -d 0 -x 20 -y 20 -z 14 -sol PW -center -l CHOL:0.446 -l PUPE:0.097    \
      -l PAPE:0.061 -l POPC:0.049 -l PUPS:0.034 -l DPPC:0.03  -l PAPS:0.028 \
      -l PAPC:0.026 -l POPE:0.025 -l POPS:0.025 -l PUPI:0.02  -l DPSM:0.015 \
      -l OUPE:0.014 -l OAPE:0.013 -l POPI:0.013 -l PAPI:0.013 -l DOPC:0.012 \
      -l PUPC:0.01  -l OUPS:0.007 -l DPPS:0.005 -l PIPI:0.005 -l DPCE:0.004 \
      -l PFPC:0.003 -l OIPC:0.003 -l OIPE:0.003 -l PNSM:0.003 -l PBSM:0.003 \
      -l PAPA:0.003 -l PAP1:0.003 -l PAP2:0.003 -l PAP3:0.003 -l PADG:0.003 \
      -l OUPC:0.002 -l POSM:0.002 -l POP1:0.002 -l POP2:0.002 -l POP3:0.002 \
      -l  IPE:0.002 -l POPA:0.001 -l DBCE:0.001 -l PNCE:0.001 -l  PPC:0.001 \
      -l  IPC:0.001 -l  PPE:0.001 -l PODG:0.001 -u CHOL:0.444 -u POPC:0.087 \
      -u DPSM:0.058 -u DPPC:0.053 -u PUPE:0.05  -u DPGS:0.049 -u PAPC:0.046 \
      -u PAPE:0.031 -u DOPC:0.022 -u PUPC:0.017 -u POPE:0.013 -u PNSM:0.013 \ 
      -u PBSM:0.011 -u PNGS:0.011 -u DPG1:0.009 -u DPG3:0.009 -u DBGS:0.009 \
      -u OAPE:0.007 -u OUPE:0.007 -u POSM:0.007 -u PFPC:0.006 -u OIPC:0.006 \
      -u POGS:0.006 -u OUPC:0.004 -u DPCE:0.004 -u PADG:0.003 -u DBG1:0.002 \
      -u PNG1:0.002 -u DBG3:0.002 -u PNG3:0.002 -u PPC:0.002  -u OIPE:0.001 \
      -u POG1:0.001 -u POG3:0.001 -u DBCE:0.001 -u PNCE:0.001 -u  IPC:0.001 \
       -u PPE:0.001 -u  IPE:0.001 -u PODG:0.001
~~~

Insane has generated a topology file for us called `GlyT2_CG_neuronal.top`.  We need to change these a little bit to add our elastic network, and include parameters for coarse grained zwitterionic glycine.

You need to edit `GlyT2_CG_neuronal.top` in a number of ways

Add `#include` statements for the files `Protein.itp` and `martini.ff/glycine-CG.itp`, and a `#define` statement to define `RUBBER_BANDS`.    Look at `GlyT2_CG.top` to see what they should look like.

Next, update the add one `Glycine` and three `NA+` molecules to the `[ molecules ]` list between the Protein and the first lipid.

If you've done it all correctly, it's time to minimise your protein in a membrane.

~~~s
gmx grompp -f emin.mdp  -c 06_GlyT2_CG_neuronal.gro  -p GlyT2_CG_neuronal.top \
           -o 07_GlyT2_neuronal_min.tpr  2>&1 | tee 07_grompp.out

gmx mdrun -v -deffnm 07_GlyT2_neuronal_min
~~~

### Neutralize and add 0.15M NaCl

insane.py has already generated solvent, but we still need to add ions, using `grompp` and `genion`.  
~~~s
gmx grompp -f emin.mdp -c 07_GlyT2_neuronal_min.gro -p GlyT2_CG_neuronal.top -o 08_GlyT2_neuronal_ions.tpr 2>&1 | tee 08_grompp.out
gmx genion -s 08_GlyT2_neuronal_ions.tpr  -p GlyT2_CG_neuronal.top -conc 0.15 -neutral -o  08_GlyT2_neuronal_ions.gro
~~~

Choose the `PW` group.

Have a look at the changes to your system and your .top file.

Gromacs has generated ions using atom and molecules names that don't quite match what Martini expects, and we need to change them from NA and CL to NA+ and CL-.  There are a lot of ions to rename, so we're going to use the program `sed` to automate the change.  To learn how `sed` works, youcan look at the man page with

~~~s
man sed
~~~

If you're a mac user, you might need to use [gsed](https://formulae.brew.sh/formula/gnu-sed) instead, or edit these commands to work with the version of `sed` that comes with Mac OS.

**Be very careful!**  `sed` is going to make a lot of changes to your files automatically, and it's not going to back them up first.  If you make a mistake, it will be very hard to undo.  Perhaps you should make a backup before you start changing things.

First, have a look at how the ions are described in `08_GlyT2_neuronal_ions.gro`.  You might need to search for `NA`

Then, run this command, and then look at how the description of the ions has changed.

~~~s
sed -i 's?NA      NA?ION    NA+?' 08_GlyT2_neuronal_ions.gro
~~~

`sed` has done a find and replace, finding `NA      NA` and replacing it with `ION    NA+`.

next, run sed three more times.  Before you run each of these commands, check what you think it's going to change, and then look to see if you were correct.

~~~s
sed -i 's?CL      CL?ION    CL-?' 08_GlyT2_neuronal_ions.gro
~~~

~~~s
sed -i 's/^NA /NA+/' GlyT2_CG_neuronal.top
~~~

~~~s
sed -i 's/^CL /CL-/' GlyT2_CG_neuronal.top
~~~

If you've done this right, it's time to minimise again!

~~~s
gmx grompp -f emin.mdp -c 08_GlyT2_neuronal_ions.gro -p GlyT2_CG_neuronal.top \
           -o 09_GlyT2_neuronal_ions_min.tpr 2>&1 | tee 09_grompp.out

gmx mdrun -v -deffnm 09_GlyT2_neuronal_ions_min
~~~

### Ready for equilibration

Once again, we need an index file.  We want to separate out system into three groups:  the `Protein`, the `lipids`, and the water and ions (`PW_ION`)

~~~s
gmx make_ndx -f 09_GlyT2_neuronal_ions_min -o GlyT2_neuronal.ndx
~~~

Look at the list of index groups.  There's a lot more than last time.

Good news:  `Protein` already exists.  So we just need to make the others.
To make `PW_ION` we're going to use an `or` statement.  Find the group numbers for the `PW` and `ION` groups, and create a new group that contains every atom that belongs to either these groups.  If ION is group 20, and PW is group 21, you would type

~~~s
20|21
~~~

and hit enter twice.

But wait, that group is called `ION_PW`.  We need to change the name to `PW_ION`. See iif you can remember how to do it with the `name` command.

To make a group of lipids, we're going to use `and` and `not` statements to find everything that isn't in the groups `Protein` or `PW_ION`

If `Protein` is group 1, and `PW_ION` is group 22, we want to type

~~~s
!1 & !22
~~~

Why do we need to use an and statement here instead of an or statement?

Oncw again, change the name of the group you have just made to `lipids`


### Equilibrate

Oncw again, we have a series of equilibration steps. You may have noticed that the membrane generated by insane is very highly ordered.

First, we're going to run equilibration with positon restraints, decreasing from 1000 to 10.  These will be different to the ones from the united atom tutorial, as we will be using a 10 femtosecond timestep.

Next, we will slowly increase the timestep from 10 femtoseconds to the final value of 20 femtosseconds.

These equilibration steps are much longer than the ones you used in your united atom tutorial, so make sure everything is ready, and then write a bash file to run the following EQ steps.

~~~s
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
~~~

Good luck!  You probably want to leave this to run overnight.

Once the EQ is complete, you should look at the trajectory of the various equilibration steps in VMD.  You'll be able to see the different lipids mixing together, and watch GlyT2 diffuse through the membrane.

### Advanced strategies: automation

This was a very long process, with a lot of steps, and to build multiple replicates you would need to repeat the entire process two more times, from start to the end of equilibration.  

Once you have a working method, where you know exactly what's going on, it's can be valuable to consider how much of your system you might want to automate.

00_README.sh is a readme file that is also a shell script.  This file contains every step from this tutorial (except downloading go_martinize.py), along with comments as to the purpose of these steps, like any good readme.  However, it is also an executable shell script.  If you run 00_README.sh, it will execute every step of this tutorial, from coarse graining GlyT2 all the way through to equilibration.

This was very valuable to me in my honours year, as I needed to use versions of this script to produce nine different simulation systems, with three replicates each.  That's twenty seven replicates total, and every one of them took at between five and eight hours to produce.  Automating this process (and thoroughly testing my automation) took me a few days, but it saved me a few weeks.

You need to be careful about deciding what to automate.  Computers are very stupid; if you don't write and test your code carefully automation provides you an opportunity to make mistakes at lighting speed and produce worthless garbage at scale.  Automating careully takes a lot of time, and if your job is small it might be better to just do it yourself.  But if you're thorough, and the job is big, sometimes automation is the right answer.  I hope considering 00_README.sh provides some insight into what is possible.