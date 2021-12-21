function R = RobustSingleRotationAveraging_from_text_file_of_x_y_z_rx_ry_rz(fn_tra_rot)
  pkg load statistics;

  R_samples = load_rotations_from_text_file_of_x_y_z_rx_ry_rz(fn_tra_rot);

  % 2-a. Average them using Hartley's L1 geodesic method 
  % (with our initialization and outlier rejection scheme):

  b_outlier_rejection = true;
  n_iterations = 10;
  thr_convergence = 0.001;
  tic;
  R_geodesic = GeodesicL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
  time_geodesic = toc;

  % 2-b. Average them using our approximate L1 chordal method 
  % (with our initialization and outlier rejection shceme)

  b_outlier_rejection = true;
  n_iterations = 10;
  thr_convergence = 0.001;
  tic;
  R_chordal = ChordalL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
  time_chordal = toc;


  % 3. Evaluate the rotation error (deg):

  %t1 = trace(R_true*R_geodesic')
  %t2 = t1 - 1
  %t3 = t2 / 2
  %t4 = acosd(t3)
  %t6 = acos(t3)
  %t5 = abs(t4)

  error_GeodesicL1Mean = abs(acosd((trace(R_true*R_geodesic')-1)/2));
  error_ChordalL1Mean = abs(acosd((trace(R_true*R_chordal')-1)/2));
  %R_geodesic
  %R_chordal
  disp(['Error (geodesic L1 mean) = ', num2str(error_GeodesicL1Mean), ' deg, took ', num2str(time_geodesic*1000), ' ms'])
  disp(['Error (chordal L1 mean) = ', num2str(error_ChordalL1Mean), ' deg, took ', num2str(time_chordal*1000), ' ms' ])


  disp('')
end



