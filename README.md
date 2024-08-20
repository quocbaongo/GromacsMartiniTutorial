# GromacsMartiniTutorial
This is a tutorial for using martini force field in Gromacs. The available material is used for coarse-grained molecular dynamics simulations workshop in the course "In silico methodologies in biochemistry and molecular medicine" provided by University of Oulu (https://opas.peppi.oulu.fi/fi/opintojakso/747613S/7534?period=2024-2025).

In this tutorial, we will use Martini force-field v3.0 (available at http://cgmartini.nl/index.php/martini-3-0) and Gromacs 2020.4 (available at https://manual.gromacs.org/documentation/2020.4/download.html). The biological system for investigation is a complex of Leptin- Leptin binding domain CRH2.

# To create environment for MDAnalysis
conda env create -f env.yml

#It is important to note that the script in the MDAnalysis directory can only be used for CGMD simulation, and it will not work for ATMD simulation. 

#The script in MDAnalysis directory requires software and python libraries specified in env.yml script.

#The backmapping approach requires python2 environment. The python, bash script and mapping files required for backmapping were obtained from http://cgmartini.nl/index.php/tutorials-general-introduction-gmx5/others-gmx5. 

