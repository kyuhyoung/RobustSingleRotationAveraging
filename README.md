# C++ version of RobustSingleRotationAveraging
OpenCV is required for matrix operations.

Checked that the outputs of c++ version are almost the same as those of original Matlab version given the same inputs. 

```console
$ make
$ ./RobustSingleRotationAveraging --th_convergence=0.001 --n_iter=10 --outlier=1
```

# Robust average of 3D transforamtion
Not only RobustSingleRotationAveraging, geometric median of 3D translation is also computed from a text file of 'x, y, z, rx, ry, rz' formatted lines.

Given the file of x_y_z_rx_ry_rz.txt each line of which is in 'x, y, z, rx, ry, rz' format,

for chordal median,
```console
$ transformation_averaging_from_text_file_of_x_y_z_rx_ry_rz x_y_z_rx_ry_rz.txt chordal
```
for geodesic median,
```console
$ transformation_averaging_from_text_file_of_x_y_z_rx_ry_rz x_y_z_rx_ry_rz.txt geodesic
```

# RobustSingleRotationAveraging
MATLAB implementation of our method proposed in "Robust Single Rotation Averaging" ([arXiv](https://arxiv.org/abs/2004.00732))

Instruction: Download our [script](https://github.com/sunghoon031/RobustSingleRotationAveraging/blob/master/RobustSingleRotationAveraging_ReleaseCode.m) and run it on Matlab.
