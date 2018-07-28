# comments
#       Y                                        _________
#       |                                       / _ \     \  
#       |                                      | / \ |     |    
#       |______Z                               | | | |     |         
#        \                                     | | | |     |                   
#         \                                    | \_/ |     |                     
#          X                                    \___/_____/         
# Left-handed axis orientation
# Positions in millimeters
# Scanner axis is z
# Use comma without space as separator in the tables. 

modality : PET
scanner name : PET_GE_DRX

scanner radius : 443   # Distance between the center of the scanner and the center of a rsector

# rsector are repeated on a ring
number of elements              : 15120
number of rsectors              : 70                         
rsectors first angle              : 90     # optional (default is 0 deg)
rsectors angular span        : 360   # optional (default is 360 deg)
rsectors ZShift                    : 0     # optional (default is none). Axial shift values must be separated by commas

number of modules transaxial    : 1     # optional (default is 1)
number of modules axial            : 4     # optional (default is 1)
module gap transaxial                : 0     # optional (default is 0mm)
module gap axial                        : 1.5 # optional (default is 0mm)

number of submodules transaxial : 1     # optional (default is 1)
number of submodules axial         : 1     # optional (default is 1)
submodule gap transaxial             : 0     # optional (default is 0mm)
submodule gap axial                     : 0     # optional (default is 0mm)

number of crystals transaxial : 9     # optional (default is 1)
number of crystals axial         : 6     # optional (default is 1)
crystal gap transaxial             : 0.065     # optional (default is 0mm)
crystal gap axial                     : 0.1      # optional (default is 0mm)

number of layers                : 1

# The 4 following parameters could defined in arrays (SizeLayer1,SizeLayer2,SizeLayer3,etc..) if the scanner includes more than one layer
crystals size depth                : 30
crystals size transaxial          : 4.23 # optional
crystals size axial                 : 6.35 # optional

# default reconstruction parameters
voxels number transaxial        : 256
voxels number axial                : 47

field of view transaxial        : 700.0
field of view axial                : 153.69

mean depth of interaction       : -1 # optional (default value : center of crystal ). Input value must correspond to the distance from the crystal surface. Negative value means default 

min angle difference        : 40 # optional (default is 0 deg)

# description
description        : Generic file of a GATE model of the GE DRX PET scanner system.
