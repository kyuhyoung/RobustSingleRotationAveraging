function rot_mat_avg = RobustSingleRotationAveraging_from_list_of_angle_axis_vectors(li_rot_vec, str_mode)
  %[li_homo_mat, li_tra, li_rot_vec] = load_translations_and_rotations_from_text_file_of_x_y_z_rx_ry_rz(fn_tra_rot);
  pkg load statistics;
  R_samples = cellfun(@RotationFromAngleAxis, li_rot_vec, 'un', 0);  
  %R_samples
  %pause(100)
  % 2-a. Average them using Hartley's L1 geodesic method 
  % (with our initialization and outlier rejection scheme):
  %str_mode
  tic;
  if (0 ~= strcmp('geodesic', lower(str_mode)))
    b_outlier_rejection = true;
    n_iterations = 10;
    thr_convergence = 0.001;
    rot_mat_avg = GeodesicL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
  elseif (0 ~= strcmp('chordal', lower(str_mode)))   
    % 2-b. Average them using our approximate L1 chordal method 
    % (with our initialization and outlier rejection shceme)
    b_outlier_rejection = true;
    n_iterations = 10;
    thr_convergence = 0.001;
    rot_mat_avg = ChordalL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
  else
    disp("Invalid input for 'str_mode'. Currently availale choices are 'geodesice' and 'chordal'");    
    return;
  endif
  time_sec = toc;  
  disp(['L1 mean of ', upper(str_mode), ' rotaion averaging took ', num2str(time_sec * 1000), ' ms'])
end