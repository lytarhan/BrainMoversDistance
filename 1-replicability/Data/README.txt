This directory contains example data used in the replicability analysis. 

Background:
In the original experiment, subjects (N = 13) viewed 120 short videos of everyday actions (more details in Tarhan & Konkle, 2019 BioRXiv; doi: https://doi.org/10.1101/618272). From this experiment, we get a beta matrix for every subject (dimensions: 120 videos x voxels).
We then reduced this matrix by selecting data only from voxels that respond reliably in the group data (using the procedure outlined in Tarhan & Konkle, 2019 BioRXiv, doi: https://doi.org/10.1101/703603).

SSBetas:
The files named (e.g.) SSBetas_Sub1.mat contain beta matrices for each subject, only in reliable voxels (reliability was calculated based on the group data to control the set of voxels across subjects).

VoxelCoordinates:
(i, j, k) coordinates for every voxel (used to make a sparsely-connected network between the voxels in order to specify a sparser constraint matrix for the Wasserstein optimization problem. This step is optional, but speeds up the process).

VoxelDistances:
square matrix of distances (in mm) between every pair of voxels. These could be measured using cortical distance, straight-line distance, or distance taking white matter connectivity into account; I chose straight-line distance.
