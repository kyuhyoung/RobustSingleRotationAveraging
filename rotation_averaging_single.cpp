#include "rotation_averaging_single.h"
#include <cmath>
#include <opencv2/sfm/numeric.hpp>

void cout_indented(int n_space, const string& str)
{
    if(n_space >= 0) std::cout << std::string(n_space * 2, ' ') << str << std::endl;
}

double lap_time(struct timespec& t_begin, struct timespec& t_end, bool is_millisecond)
{
    double lap_sec = (t_end.tv_sec - t_begin.tv_sec) + (t_end.tv_nsec - t_begin.tv_nsec) / 1000000000.0;
    return is_millisecond ? 1000.0 * lap_sec : lap_sec;
}

float acosd(float cosain)
{
    //float acos_rad = acos(cosain);
    return 180.0 * acos(cosain) / CV_PI;     
}

Mat zeros_like(const Mat& mat_ori)
{
    return Mat::zeros(mat_ori.size(), mat_ori.type());
}

template <typename T> 
T median_of_vector(vector<T> &v, int n_sp)
{
    cout_indented(n_sp, "median_of_vector START");
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin()+n, v.end());
    cout_indented(n_sp + 1, "median : " + to_string(v[n]));
    cout_indented(n_sp, "median_of_vector END");
    return v[n];
}

Mat ProjectOntoSO3(const Mat& rot_mat_ori)
{
    return rot_mat_ori;
}    

Mat logarithm_map(const Mat& rot_mat)
{
    Mat rot_vec;
    Rodrigues(rot_mat, rot_vec);
    return rot_vec;    
}


Mat SkewSymmetricMatrix(const Mat& vec)
{
    return cv::sfm::skew(vec);
} 

Mat RotationFromUnitAxisAngle(const Mat& unit_axis, double angle)
{
    Mat rot_mat;
    Rodrigues(angle * unit_axis, rot_mat);
    return rot_mat;
}    

Mat GeodesicL1Mean(vector<Mat> R_input, bool b_outlier_rejection, int n_iterations, float thr_convergence)
{
    int n_samples = R_input.size();
    int idx_firstQ = ceil(float(n_samples) / 4.0);
    Mat s = zeros_like(R_input[0]);
    for(int iR = 0; iR < 3; iR++)
    {
        for(int iC = 0; iC < 3; iC++)
        {
            vector<double> li_val(n_samples);
            for(int iS = 0; iS < n_samples; iS++)
            {
                li_val[iS] = R_input[iS].at<double>(iR, iC);
            }
            double val_med = median_of_vector<double>(li_val, -100);
            s.at<double>(iR, iC) = val_med;
        }
    }     
    
    Mat R = ProjectOntoSO3(s);
    for(int iI = 0; iI < n_iterations; iI++)
    {
        vector<Mat> vs(n_samples);
        vector<double> v_norms(n_samples);
        for(int iS = 0; iS < n_samples; iS++)
        {
            Mat v = logarithm_map(R_input[iS] * R.t());
            double v_norm = norm(v);
            vs[iS] = v;
            v_norms[iS] = v_norm;
        }
        double thr = 100000000000000000000.0;
        if(b_outlier_rejection)
        {
            vector<double> sorted_v_norms(v_norms);
            sort(sorted_v_norms.begin(), sorted_v_norms.end());
            double v_norm_firstQ = sorted_v_norms[idx_firstQ];
            thr = n_samples <= 50 ? MAX(v_norm_firstQ, 1.0) : MAX(v_norm_firstQ, 0.5);
        }
        Mat step_num = zeros_like(vs[0]);
        double step_den = 0;
        for(int iS = 0; iS < n_samples; iS++)
        {
            Mat v = vs[iS];
            double v_norm = v_norms[iS];
            if(v_norm > thr)
            {
                continue;
            }
            step_num += v / v_norm;
            step_den += 1.0 / v_norm;
        }
        Mat delta = step_num / step_den;
        double delta_angle_rad = norm(delta);
        Mat delta_axis = delta / delta_angle_rad;
        Mat R_delta = RotationFromUnitAxisAngle(delta_axis, delta_angle_rad);
        R = R_delta * R;
        if(delta_angle_rad < thr_convergence)
        {
            break;
        }
    }
    return R;
}

Mat avoid_median_to_be_the_same_as_one_of_samples(const vector<Mat>& R_input, const Mat& s_ori)
{
    bool is_current_median_the_same_as_one_or_more_of_samples = false;
    for(auto R : R_input)
    {
        double sum_dif = norm(R - s_ori);
        if(sum_dif < 0.001)
        {
            is_current_median_the_same_as_one_or_more_of_samples = true;
            break;
        }     
    }
    if(is_current_median_the_same_as_one_or_more_of_samples)
    {
        Mat mat_rand = zeros_like(s_ori);
        randu(mat_rand, Scalar(-0.001), Scalar(0.001));
        Mat s_new = s_ori + mat_rand;
        return s_new;
    }
    else
    {
        return s_ori;
    }
}

Mat ChordalL1Mean(vector<Mat> R_input, bool b_outlier_rejection, int n_iterations, float thr_convergence)
{
    int n_samples = R_input.size();
    int idx_firstQ = ceil(float(n_samples) / 4.0);
    Mat s = zeros_like(R_input[0]);
    for(int iR = 0; iR < 3; iR++)
    {
        for(int iC = 0; iC < 3; iC++)
        {
            vector<double> li_val(n_samples);
            for(int iS = 0; iS < n_samples; iS++)
            {
                li_val[iS] = R_input[iS].at<double>(iR, iC);
            }
            double val_med = median_of_vector<double>(li_val, -100);
            s.at<double>(iR, iC) = val_med;
        }
    }     
 
    for(int iI = 0; iI < n_iterations; iI++)
    {
        s = avoid_median_to_be_the_same_as_one_of_samples(R_input, s);
        vector<double> v_norms(n_samples);
        for(int iS = 0; iS < n_samples; iS++)
        {
            Mat v = R_input[iS] - s;
            double v_norm = norm(v);
            v_norms[iS] = v_norm;
        }
        
        double thr = 100000000000000000000.0;
        if(b_outlier_rejection)
        {
            vector<double> sorted_v_norms(v_norms);
            sort(sorted_v_norms.begin(), sorted_v_norms.end());
            double v_norm_firstQ = sorted_v_norms[idx_firstQ];
            thr = n_samples <= 50 ? MAX(v_norm_firstQ, 1.356) : MAX(v_norm_firstQ, 0.7);
            //  2*sqrt(2)*sin(1/2) is approximately 1.356
            //  2*sqrt(2)*sin(0.5/2) is approximately 0.7
        }
        Mat step_num = zeros_like(s);
        double step_den = 0;
        for(int iS = 0; iS < n_samples; iS++)
        {
            double v_norm = v_norms[iS];
            if(v_norm > thr)
            {
                continue;
            }
            step_num += R_input[iS] / v_norm;
            step_den += 1.0 / v_norm;
        }
        Mat s_prev = s;
        s = step_num / step_den;
        Mat update_medvec = s - s_prev;
        if(norm(update_medvec) < thr_convergence)
        {
            break;
        }     
    }     
    Mat R = ProjectOntoSO3(s);
    return R;
}

pair<Mat, vector<Mat> > prepare_noisy_rotation_matrix()
{
    Mat R_true = (Mat_<double>(3, 3) << 0.880924, 0.277119, 0.383637, 0.013324, 0.795787, -0.605429, -0.473070, 0.538449, 0.697336);
    vector<Mat> R_samples(7);
    R_samples[0] = (Mat_<double>(3, 3) << 0.87568, 0.30871, 0.37134, -0.015865, 0.78695, -0.61681, -0.48264, 0.53423, 0.69401);           
    R_samples[1] = (Mat_<double>(3, 3) << 0.88084, 0.27657, 0.38424, 0.014008, 0.79603, -0.60509, -0.47321, 0.53837, 0.6973);
    R_samples[2] = (Mat_<double>(3, 3) << 0.89546, 0.2712, 0.35299, 0.012714, 0.77708, -0.62928, -0.44496, 0.56798, 0.69239);
    R_samples[3] = (Mat_<double>(3, 3) << 0.85608, 0.27909, 0.43502, 0.060279, 0.78202, -0.62033, -0.51332, 0.55727, 0.65265);
    R_samples[4] = (Mat_<double>(3, 3) << 0.87326, 0.26435, 0.40932, 0.030012, 0.80926, -0.58668, -0.48634, 0.52461, 0.69876);
    R_samples[5] = (Mat_<double>(3, 3) << -0.33503, 0.24321, -0.91028, -0.061857, 0.95835, 0.27882, 0.94017, 0.14972, -0.30603);
    R_samples[6] = (Mat_<double>(3, 3) << -0.073099, 0.96268, 0.26058, 0.9324, 0.1587, -0.32472, -0.35396, 0.21923, -0.9092);
    return make_pair(R_true, R_samples);
}


