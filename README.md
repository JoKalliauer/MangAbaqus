# MangAbaqus

Prof. Mang's Project using Abaqus-FEM, originally created by Michał Malendowski (©2019-2020) and extended and rebuilt by Johannes Kalliauer (©2020-2023).

This code was published in Aug. 2022 under GNU AGPLv3 for the paper _Conditions for minimum stiffness of proportionally loaded structures_ .

To start call Main_*.m
- eccentric Compression:  Main_ecc.m
- Bending: Main_bend.m
- Stützlinienbogen: Mail_TLArch.m
- Cantilever: Main_canti.m

This code was used the following publiations
- [On a remarkable geometric-mechanical synergism based on a novel linear eigenvalue problem](https://doi.org/10.1007/s00707-021-03091-5)
- [Are the terms stiffening/softening structures mechanically unambiguous?](https://doi.org/10.1016/j.euromechsol.2022.104756)
- [Conditions for minimum stiffness of proportionally loaded structures](https://doi.org/10.1016/j.cma.2022.115820)

This projects were financial supported by the Austrian Science Fund (FWF) in the framework of the research project P 31617-N32 [Pseudo-kinematic invariants - gems in FE structural analyses].


## Basic Structure

Several Files to start Matlab:  
- Main_TLArch.m ... for Trust-line-arch
- Main_bend ... for pure bending of a beam on two supports
- Main_ecc_ACME ... for excentric compression according to Acta Mechanica-Paper (bending about the strong axis)
- Main_ecc ... excentric compression of the newest version (bending about the weak axis)
- ...

Those Main_*.m files call Abaqus_single_run.m which calls the following functions (simplified)
* runEigenProblem ... run the Eingenvalue-Problem
  * selectModel .. calls a function to create the input-file 
  * AbaqusModelsGeneration.runAbaqus ... run the input-file in Abaqus
  * AbaqusModelsGeneration.getStiffnessMatrices ... get the stiffness-matrix from Abaqus-results
  * AbaqusModelsGeneration.getHistoryOutputFromDatFile ... get the nodal-results from Abaus
  * runEigenProblemSub ... Run the core of the eigenvalue-Problem and saving it into model
    * solveCLEforMinEigNew ... get one specific eigenvalue and the eigenvector
  * runEigenProblemDispJK ... Posprocessing the displacements
* sortEigenValuesAndGetQuantities ... does the calculation of \rho
* plotresMulti ... Plots the requested graphs
 

## Matlab-Files

- Abaqus_single_run ... the function to call, which calls the individual subprograms
- compareEigenmodes
- eccfromU ... calculate the exccentricity of a compressed beam for a given energy-ratio
- getQuantities ... calulate rho, acceleration, ...
- getQuantities3D ... Convert N-dimensional to 3-Dimensional and then calculate rho (outdated)
- getQuantitiesStern ... calulate rho based on rho (outdated)
- GetSize ... 
- InterpolateJK ... does a linear interploation
- JKExtend
- Main_Abaqus_single_run ... oldest Version for running any example
- Main_bend ... parameters for pure bending
- Main_bock ... snap-through of two beams, load-controlled
- Main_bockDisp  ... snap-through of two beams, displacement-controlled
- Main_canti .. parameters for cantilever beam
- Main_detKt ... example for Fig6 in https://doi.org/10.1016/j.euromechsol.2022.104756
- Main_ecc ... exzentric compression, weak axis
- Main_ecc_ACME ... excentric compression, strong axis
- Main_TLArch ... prarameters for trust-line-arch
- markMins ... mark miniums in a graph and determine a usefull upper limit of the y-axis of the graph
- markWende ... mark infelction points in a graph
- mergeModel ... if loadfactor=0, then the model runs it with loadfactor=-1 and +1, and then merges both to one
- myIsField ... checks if fieldname is part of structure
- NodalResults2Displ  ... Returns the displacements (inkl. rotations) based on the nodal results
- NodalResults2DisplJK ... Splits the Abaqus-results in displacements and rotations
- PlotCSV2 ... plot graph for paper
- PlotCSV2TL ... plot graph for trust-line-arch
- PlotCSVPaper  ... plot graph for paper
- PlotCSVPaperIso ... plot isometric-view
- PlotCSVPaperold ... deleted
- plotitJK ... plot with predefined settings by Johannes Kalliauer
- plotres ... plot the results (old by Malendowski)
- plotresMulti ... plot several results based on plotfig
- printresMulti ... save results in a spreadsheat
- Profil ... calulate profile-dimensions
- RealEV ... split the Eivenvetor in a matrix (dof-per-node x nodes)
- runEigenProblem ... the program that calculates the eigenvalues
- runEigenProblemDispJK ... subprogram for deriving the displacements
- runEigenProblemRotJK ... subprogram for deriving the rotations
- runEigenProblemSub ... the core suprogram of runEigenProblem
- selectModel ... call the right function based on 'modelprops.testcase'
- solveCLEforMinEigNew ... solve the eigenvalue-problem
- sortEigenValuesAndGetQuantities ... normizaiton of the vector for deriving rho


## Git-Files
*	.gitattributes … see https://dev.to/deadlybyte/please-add-gitattributes-to-your-git-repository-1jld
*	.gitignore … which file-types should not be saved on git
*	CHANGELOG ... recent changes to the project
*	LICENSE ... tells you how to reuse the project
*	README.md … This explantation
*	gitpush.sh … A code for git-pushing from Linux
*	gitpush.bat … A code for git-pushing from Windows