modality: SPECT_CONVERGENT
scanner name: SPECT_CASTOR_BENCHMARK
description: CASToR SPECT benchmark scanner for versions 1.X

number of detector heads: 2

trans number of pixels: 1 # 1 in case of monolythic
trans pixel size: 613.7856 # Given in mm
trans gap size: 0 # Given in mm

axial number of pixels: 1
axial pixel size: 613.7856 # Given in mm
axial gap size: 0 # Given in mm

detector depth: 20

# Distance between the center of rotation (COR) of the scanner and the surface of a detection head in mm
scanner radius: 142.60, 177.10 # Head 1 and 2

# Collimator configuration

head1:
trans focal model: constant
trans number of coef model: 1
trans parameters: 0 # focal distance in mm, 0 for parallel
axial focal model: constant
axial number of coef model: 1
axial parameters: 0

head2:
trans focal model: constant
trans number of coef model: 1
trans parameters: 0 # focal distance in mm, 0 for parallel
axial focal model: constant
axial number of coef model: 1
axial parameters: 0

voxels number transaxial: 128 # optional (default is the half of the scanner radius)
voxels number axial: 128 # optional (default is length of the scanner computed from the given parameters)

field of view transaxial: 613.7856 # optional (default is the half of the scanner radius)
field of view axial: 613.7856 # optional (default is length of the scanner computed from the given parameters)
