# GromacsMartiniTutorial
This is a tutorial for using martini force field in Gromacs. The available material is used for coarse-grained molecular dynamics simulations workshop in the course "In silico methodologies in biochemistry and molecular medicine" provided by University of Oulu (https://opas.peppi.oulu.fi/fi/opintojakso/747613S/7534?period=2024-2025).

In this tutorial, we will use Martini force-field v3.0 (available at http://cgmartini.nl/index.php/martini-3-0) and Gromacs 2020.4 (available at https://manual.gromacs.org/documentation/2020.4/download.html). The biological system for investigation is a complex of Leptin- Leptin binding domain CRH2.

## 1. Simulation environment
To create the required environment for performing coarse-grained molecular dynamics (CGMD) simulations and their subsequent analyses, use the following command: conda env create -f env.yml

It is important to note that the script in the MDAnalysis directory can only be used for CGMD simulation, and it will not work for ATMD simulation. 

The script in MDAnalysis directory requires software and python libraries specified in env.yml script.

The backmapping approach requires python2 environment. The python, bash script and mapping files required for backmapping were obtained from http://cgmartini.nl/index.php/tutorials-general-introduction-gmx5/others-gmx5. 

## 2. Tutorial structure
Detailed instructions how to execute a coarse-grained molecular dynamics simulations (CGMD) using Martini v3.0 force-field can be found in the file 'MartiniCG_Execution_Analysis.pptx'

Directory 'mdp' contains all the .mdp files that are required for this tutorial. All the files were downloaded and modified where it is necessary from the tutorial provided by Martini developers (http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/564-2-proteins-basic-and-martinize-2)

Directory 'CG_with_elasticNetwork' and 'CG_with_NoElasticNetwork' contains the pre-generated topologies and trajectories of the Leptin-LR CRH2 complex simulations with and without the presence of elastic network, respectively.

Directory 'MDAnalysis' contains the python scripts used for trajectory analysis. Invoking 'python3 file_name --help' to access the description of the program and instruction, how to use the program.

Directory 'BackMappedCGtoAT' contains all the required material for transforming the protein structure in coarse-grained resolution to all-atom resolution (check http://www.cgmartini.nl/index.php/tutorials/37-tutorial2/314-tutorial-reverse-mapping for further information).

Directory 'Student_feedback_2024' includes detailed student feedback for the course in 2024
