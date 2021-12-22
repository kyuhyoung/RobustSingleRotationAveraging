# C++ version of RobustSingleRotationAveraging
OpenCV is required for matrix operations.

Checked that the outputs of c++ version are almost the same as those of original Matlab version given the same inputs. 

```console
$ make
$ ./RobustSingleRotationAveraging --th_convergence=0.001 --n_iter=10 --outlier=1
```

# RobustSingleRotationAveraging
MATLAB implementation of our method proposed in "Robust Single Rotation Averaging" ([arXiv](https://arxiv.org/abs/2004.00732))

Instruction: Download our [script](https://github.com/sunghoon031/RobustSingleRotationAveraging/blob/master/RobustSingleRotationAveraging_ReleaseCode.m) and run it on Matlab.
