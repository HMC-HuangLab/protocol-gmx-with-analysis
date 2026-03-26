# **Automated MD Analysis Pipeline for Protein-Ligand Systems**

This repository provides a high-performance, automated pipeline for analyzing Protein-Ligand Molecular Dynamics (MD) trajectories. It supports multi-replica processing, Free Energy Surface (FES) construction, and publication-quality (**1000 DPI**) image generation.

## **📜 Acknowledgments**

This analysis suite is built upon the foundational work of **Thibault Tubiana, PhD**, specifically his repository: [**protocolGromacs**](https://github.com/tubiana/protocolGromacs) (*My personal GROMACS workflow and best practices*). We sincerely thank him for his contribution to the GROMACS community and for providing a robust protocol architecture.

## **🚀 1\. GMX Environment Setup**

It is highly recommended to use conda to create an isolated environment to ensure compatibility between gmx\_MMPBSA, AmberTools, and GROMACS.

\# Create and activate environment  
conda create \-n gmx python=3.9  
conda activate gmx

\# Install core dependencies (Topology & Chemistry)  
conda install \-c conda-forge \-c salilab acpype  
conda install \-c conda-forge openbabel  
conda install pdbfixer  
conda install ipykernel  
conda install biopython jupyter nbconvert 

\# Install calculation tools  
conda install \-c ostrokach dssp  
conda install \-c conda-forge mpi4py=3.1.3 ambertools=21.12   
conda install cuda-cudart cuda-version=12

\# Install specific ParmEd version (Critical for gmx\_MMPBSA)  
python \-m pip install git+\[http://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4\](http://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4)

\# Install GROMACS (2021.3 with GPU support)  
conda install \-c bioconda gromacs=2021.3 

\# Install gmx\_MMPBSA  
python \-m pip install gmx\_MMPBSA

\# Install Visualization tools  
conda install pymol-open-source

## **🏗️ 2\. Pipeline Architecture**

The analysis script anal\_test.sh follows a modular workflow to handle multiple simulation replicas:

1. **Trajectory Standardization**: Synchronizes topologies and cleans trajectories (water removal, complex centering) to generate match.xtc and match.tpr.  
2. **Global Metrics**: Calculates RMSD (Protein & Ligand), RMSF, Radius of Gyration (Rg), and SASA.  
3. **Replica-Specific Metrics**: Individual Hydrogen Bond evolution and 2D Free Energy Surfaces (FES) with colorbar scales.  
4. **Advanced Thermodynamics**: Binding Free Energy calculations via gmx\_MMPBSA with Interaction Entropy (IE) corrections.  
5. **Residue-Level Deep Mining**:  
   * Automatically parses all residues (34+) contributing to the Binding Energy (DELTAS).  
   * Calculates Residue-Ligand **Min-Distance** and **Occupancy**.  
6. **GUI Replication Factory**:  
   * Automatically generates 5 core TIF images (TDC, SDC, GB Delta, etc.) replicating the gmx\_MMPBSA\_ana GUI style.  
   * **High Resolution**: All images are exported at **1000 DPI** for publication.

## **🧪 3\. Default Force Field & Parameters**

The pipeline adheres to the following scientific standards:

* **Protein**: Amberff14SB (Optimized for kinase and complex systems).  
* **Ligand**: GAFF2 (Generated via acpype with BCC charges).  
* **Water Model**: TIP3P.  
* **MMPBSA Settings**:  
  * igb \= 2 / 5 (Generalized Born models).  
  * idecomp \= 2 (Per-residue energy decomposition).  
  * saltcon \= 0.15 (Physiological salt concentration).

## **🛠️ 4\. Usage**

### **Execution**

Use numactl to bind specific CPU cores for maximum efficiency during parallel MMPBSA tasks:

\# Syntax: ./anal\_test.sh \[Workflow\_Path\] \[Ligand\_Name\] \[GPU\_ID\]  
numactl \--physcpubind=52-67 ./anal\_test.sh /path/to/workdir LIG 2

### **Folder Structure (Results Archive)**

All results are automatically organized within the merged\_analysis directory:

merged\_analysis/  
├── rmsd\_merged.csv            \# Comparison data across all replicas  
├── replica\_1/                 \# Specific data for Replica 1  
│   ├── distance\_plots/        \# 1000 DPI plots for 34+ individual residues  
│   ├── hbond\_rep1.png         \# Hydrogen bond count over time  
│   ├── fes2d\_rep1.png         \# 2D Free Energy Surface with Colorbar  
│   ├── summary\_occupancy.png  \# Bar chart of contact frequencies  
│   └── cpd3\_rep1\_SDC.tif      \# GUI-style Sidechain Decomposition (TIF)  
└── replica\_2/                 \# Data specific to Replica 2 ...

## **📌 5\. Notes**

* **DPI Settings**: All script-generated figures are forced to **1000 DPI**. TIF files may be large; ensure adequate disk space.  
* **Naming Convention**: The script automatically handles GLY:193 colon-separated residue naming from MMPBSA output.

# **Original protocolGromacs Documentation**

*Author: Thibault Tubiana, PhD* **Please read before using this script.**

## **Description**

This script is made to facilitate the preparation and production of protein and protein/ligand, MD.

It follows the procedure described for teaching I made at the University of Bergen. You can find lectures content on this page http://tubiana.me/teaching/kjem220-molecular-modelling/ or the pdf describing all the steps on this script here: http://tubiana.me/teaching\_files/biocat2020/Tutorial\_Gromacs-2019.pdf

Fundamental analysis is also generated with gromacs tools (temperature/pressure/rmsd/rmsf/...), and the production trajectories are also cleaned with trjconv (imaging/protein centred/water stripped), but all the original trajectories are kept.

Feel free to make other analysis of course, like trajectory clustering with TTClust https://github.com/tubiana/TTClust 😇

## **Disclamer**

* Each system is unique. This protocol and MD parameters is not adapted for all systems. If your system crash, you may have to tweak MDP parameters.  
* Ligand parametrisation is *"quick and dirty"*, For a more stable MD system you may have to tweak ACPYPE parameters (and check the hydrogens that are added with babel).

## **How to**

1. Make sure you have all the dependencies  
   1. If you have a protein-ligand system, make sure acpype is installed  
   2. Gromacs  
   3. (optional) DSSP version 3  
2. Clone this repository with the command git clone https://github.com/tubiana/protocolGromacs.git  
3. Put your PDB in the repository  
4. Make the change you need in runGromacs.sh  
5. run the script with bash runGromacs.sh

## **Parameters**

* **FILE**: PDB filename without the extension  
* **LIGNAME**: 3 letter ligand name  
* **NUMBEROFREPLICAS**: Number of replicas (Default is 3\)  
* **SIMULATIONTIME**: Simulation time in ns.

## **Last thing...**

Have fun with MD and send me a mail if or open an issue if you have any problems, or just if you used this script and want to thanks me, I will be please to know that it was useful for someone 🙂

Thibault Tubiana.
