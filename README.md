# pose_trees
This repository includes implementations that were development and evaluated for a Bachelors Thesis, which is included above(soon). The implementation where build upon an Octree/Quadtree/N-dimensional linear tree library from [attcs/Octree](https://github.com/attcs/Octree).
All inplementation are inlcuded as header files for a 3D Pose. This repo also includes the file that was used to create the data of the results. This project was compiled with the GCC compiler and the CmakeList should work for that. There is also a package.xml included in case one want to work with this project under ROS2, the CMake propably needs some modifications then.

## Short Overview see Thesis for more information
This repository includes implementations to test NNS times in N-dimensional axis-alligned trees for poses(position + orientation). This includes 5 approaches to with different aspects to compensate the continouty property of Euler Angles.

Add Image...

Best approaches...
Use quaternion or AlgoV2 depending on the usecase.
