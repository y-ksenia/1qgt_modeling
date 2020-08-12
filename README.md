# Investigation of the mechanical properties of viruses in silico

This is a repository with the code of my BSs (2017) and MSc (2019) theses at the Moscow Institute of Physics and Technology 
in Laboratory for Computer and Mathematical Modeling of Biological Systems "In silico investigation of the mechanical properties of DNA-containing capsids". 
You can find the full final text of my research as MSc_thesis_MIPT.pdf.

## Structure
#### Creation of structures and their topologies
- build_capsid - Creation of the Hepatitis B capsid (https://www.rcsb.org/pdb/explore/remediatedSequence.do?structureId=1QGT)
- create_topology_dna - DNA structure and topology creation scripts
- change_ion - scripts for changing the concentration of ions in the box
- parse_dcd - additional scripts to eliminate errors in dsd-files of particle trajectories
#### Analysis
- calc_pers_length - calculation of persistence length of DNA
- calc_dna_central_line - calculating the position of the central line of double-stranded DNA
- calc_interchain_dist - calculation of DNA interchain distances using central line
- calc_rad_distr - calculation of radial distribution of DNA packed into the capsid
#### Visualization
- build_pressure_plots - visualization of mechanic properties for the experiments
