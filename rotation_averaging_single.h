#ifndef _ROTATION_AVERAGING_SINGLE_H_
#define _ROTATION_AVERAGING_SINGLE_H_

#include <opencv2/opencv.hpp>
#include <sys/time.h>
using namespace cv;
using namespace std;

double lap_time(struct timespec& t_begin, struct timespec& t_end, bool is_millisecond);

double acosd(double cosain);

Mat GeodesicL1Mean(vector<Mat> li_rot_mat, bool outlier_rejection, int n_iterations, float thr_convergence);
Mat ChordalL1Mean(vector<Mat> li_rot_mat, bool outlier_rejection, int n_iterations, float thr_convergence);
    
Mat zeros_like(const Mat& mat_ori);

Mat logarithm_map(const Mat& rot_mat);

Mat SkewSymmetricMatrix(const Mat& vec); 

Mat RotationFromUnitAxisAngle(const Mat& unit_axis, double angle);

pair<Mat, vector<Mat> > prepare_noisy_rotation_matrix();

#endif //   _ROTATION_AVERAGING_SINGLE_H_
