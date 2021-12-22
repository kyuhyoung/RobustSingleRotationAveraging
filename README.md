# C++ version of RobustSingleRotationAveraging
OpenCV is required for matrix operations.

Checked that the outputs of c++ version are almost the same as those of original Matlab version given the same inputs. 

```console
$ make
$ ./RobustSingleRotationAveraging --th_convergence=0.001 --n_iter=10 --outlier=1
```

# Robust averaging of 3D transforamtion (Matlab)
Not only RobustSingleRotationAveraging, geometric median of 3D translation is also computed from a text file of 'x, y, z, rx, ry, rz' formatted lines.

Given the file 'x_y_z_rx_ry_rz.txt' each line of which is in 'x, y, z, rx, ry, rz' format,

in Matlab Command Window, for chordal median,
```console
>> transformation_averaging_from_text_file_of_x_y_z_rx_ry_rz x_y_z_rx_ry_rz.txt chordal
```
in Matlab Command Window, for geodesic median,
```console
>> transformation_averaging_from_text_file_of_x_y_z_rx_ry_rz x_y_z_rx_ry_rz.txt geodesic
```

![image](https://user-images.githubusercontent.com/12492992/147081735-52fc8f1a-57bc-443c-bd92-1e02edbf738c.png)

![image](https://user-images.githubusercontent.com/12492992/147081816-50f65f62-9c03-4607-9261-b9c7d5f4ae57.png)


# RobustSingleRotationAveraging
MATLAB implementation of our method proposed in "Robust Single Rotation Averaging" ([arXiv](https://arxiv.org/abs/2004.00732))

Instruction: Download our [script](https://github.com/sunghoon031/RobustSingleRotationAveraging/blob/master/RobustSingleRotationAveraging_ReleaseCode.m) and run it on Matlab.
