#include "rotation_averaging_single.h"

cv::String keys =
    //"{@image            |<none>             | input image path}"         // input image is the first argument (positional)
    //"{@face             |/path/to/file.xml  | face cascade path}"        // optional, face cascade is the second argument (positional)
    "{idx_init          |0                  | Initial index of image file or pose file}"
    "{n_iter            |                  | Number of iteration in optimization}"
    "{th_converge       |              | Small enough angle (radian) difference to consider to be converged in optimization step.}"
    "{is_robot_eye      |0                  | Zero if hand-eye calibration. Non-zero if robot-eye calibration.}"         // optional, default value ""
    "{input_mode        |                   | Should be one of 'pnp', 'image', 'image+p2d' and 'p2d'.}"
    "{dir_p2d           |                   | Directory of Matlab toolbox corner det results}"
    "{p2d_refine_he     |                  | Use Matlab toolbox corner det results for seed of chessboard refinement in HE if non-zero}"
    "{is_verbose        |0                  | Print out optimization process during st-handeye interation}"
    "{outlier        |1                  | Reject outlier if this is non-zero. Otherwise do not apply outlier rejection}"
    "{robot_type        |UR                | Robot type}"         // optional, default value ""
    "{img_ext           |bmp                | Image file extension}"         // optional, default value ""
    "{prefix_corner     |                   | Prefix of Matlab corner files}"
    "{yml_cam           |                   | Path to camera instrinsic parameter yml file}"
    "{is_calib_test     |                  | Zero if it is for handeye or roboteye calibration. Nonzero if it is for calibration test}"
    "{path_img_chessboard    |                   | Path to chessboard image }"
    "{mm_gripper            |                   | Physical length of gripper in millimeter}"
    "{mm_tip            |                   | Physical length of additional tip mounted on gripper in millimeter}"         // optional, default value ""
    "{help      |      | show help message}";      // optional, show help optional
int main( int argc, char* argv[] )
{
#if 0
    double angle_perturb = 0.032111;
    Mat axis_perturb = (Mat_<double>(3, 1) << -0.38128, 0.46157, 0.80099);    
    Mat angle_axis_purburb = angle_perturb * axis_perturb;
    Mat R_perturb;
    Rodrigues(angle_axis_purburb, R_perturb);
    Mat R_perturb_2 = RotationFromUnitAxisAngle(axis_perturb, angle_perturb);
    cout << "angle_perturb : " << angle_perturb << endl;
    cout << "axis_perturb : " << axis_perturb.t() << endl;
    cout << "angle_axis_purburb : " << angle_axis_purburb.t() << endl;
    cout << "R_perturb : " << endl << R_perturb << endl;
    cout << "R_perturb_2 : " << endl << R_perturb_2 << endl;
    Mat rot_mat = (Mat_<double>(3, 3) << 0.9995373, 0.0281356, 0.0120115, -0.0281770, 0.9995956, 0.0033443, -0.0119148, -0.0036842, 0.9999187);

    Mat rot_vec = logarithm_map(rot_mat);
    cout << "rot_mat : " << endl << rot_mat << endl;
    cout << "rot_vec : " << rot_vec.t() << endl;
    Mat vec = (Mat_<double>(3, 1) << 0.78711, 0.58950, -0.18152);
    Mat mat = SkewSymmetricMatrix(vec);
    cout << "vec : " << vec.t() << endl;
    cout << "mat : " << endl << mat << endl;
    exit(0);
#endif  //  1    
    
    CommandLineParser parser( argc, argv, keys );
    parser.about( "rotation averaging" );
    if ( parser.has( "help" ) ) { parser.printMessage();  return 0; }
#if 0    
    bool is_robot_eye = 0 != parser.get<int>( "is_robot_eye" );
    std::string path_calib_yml = parser.get<std::string>( "path_calib_yml" );
    bool is_calib_test = 0 != parser.get<int>( "is_calib_test" );
    std::string yml_cam = parser.get<std::string>( "yml_cam" );
    std::string yml_chessboard = parser.get<std::string>( "yml_chessboard" );
    std::string yml_bMc_or_eMc = parser.get<std::string>( "yml_bMc_or_eMc" );
    std::string yml_bMe = parser.get<std::string>( "yml_bMe" );
    std::string yml_eTt_meter = parser.get<std::string>( "yml_eTt_meter" );
    std::string path_img_chessboard = parser.get<std::string>( "path_img_chessboard" );
    float mm_gripper = parser.get<float>( "mm_gripper" ), mm_tip = parser.get<int>( "mm_tip" );
    std::string str_robot_type = parser.get<std::string>( "robot_type" );

#endif  //  0    
    float thr_convergence = parser.get<float>( "th_converge" );
    int n_iterations = parser.get<int>( "n_iter" );
    bool b_outlier_rejection = 0 != parser.get<int>( "outlier" );
    pair<Mat, vector<Mat> > pa_R_true_R_samples = prepare_noisy_rotation_matrix();
    Mat R_true = pa_R_true_R_samples.first;
    vector<Mat> R_samples = pa_R_true_R_samples.second;
    struct timespec t_begin, t_end_1, t_end_2;
    clock_gettime( CLOCK_MONOTONIC, &t_begin);
    Mat R_geodesic = GeodesicL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
    clock_gettime( CLOCK_MONOTONIC, &t_end_1);
    #if 0
    Scalar t1 = trace(R_true * R_geodesic.t());
    cout << "t1 : " << t1 << endl;
    double t2 = t1.val[0] - 1.0;
    cout << "t2 : " << t2 << endl;
    double t3 = t2 / 2.0;
    cout << "t3 : " << t3 << endl;
    double t4 = acosd(t3);
    cout << "t4 : " << t4 << endl;
    double t5 = fabs(t4);
    cout << "t5 : " << t5 << endl;
    #endif  //  0
    float error_GeodesicL1Mean = fabs(acosd((trace(R_true * R_geodesic.t()).val[0] - 1.0) / 2.0));
    double time_geodesic = lap_time(t_begin, t_end_1, true);
    cout << "Error (geodesic L1 mean) = " << error_GeodesicL1Mean << " deg, took " << time_geodesic << " ms" << endl;
    //exit(0);
    Mat R_chordal = ChordalL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
    clock_gettime( CLOCK_MONOTONIC, &t_end_2);
    //cout << "R_geodesic : " << endl << R_geodesic << endl;    
    //cout << "R_chordal : " << endl << R_chordal << endl;    
    //error_GeodesicL1Mean = abs(acosd((trace(R_true*R_geodesic')-1)/2));
    //error_ChordalL1Mean = abs(acosd((trace(R_true*R_chordal')-1)/2));
    float error_ChordalL1Mean = fabs(acosd((trace(R_true * R_chordal.t()).val[0] - 1.0) / 2.0));
    double time_chordal = lap_time(t_end_1, t_end_2, true);
    cout << "Error (chordal L1 mean) = " << error_ChordalL1Mean << " deg, took " << time_chordal << " ms" << endl;
    return 1;
}  
