function tra_avg = median_translation_from_list_of_translation_vectors(li_tra_cell)
  %[li_homo_mat, li_tra_cell, li_rot_vec] = load_translations_and_rotations_from_text_file_of_x_y_z_rx_ry_rz(fn_tra_rot);   

  li_tra_mat = [];
  for iC = 1 : length(li_tra_cell)
    li_tra_mat(:, end + 1) = li_tra_cell{iC};
  endfor
  %class(li_tra_mat)
  %pause(100);
  tic;
  if 0
    tra_med = geometric_median(li_tra_mat);
    time_med = toc;
    %tra_med'
    disp(['Geometric median of translation took ', num2str(time_med * 1000), ' ms'])
    tra_avg = tra_med;
  else
    %[mu, S] = DetectMultVarOutliers(li_tra_mat', [], [], false);
    [tra_outlier, S, RD, chi_crt] = DetectMultVarOutliers(li_tra_mat', [], [], false);
    time_outlier = toc;
    %tra_outlier
    %S
    %RD
    %chi_crt
    disp(['Outlier aware median of translation took ', num2str(time_outlier * 1000), ' ms'])  
    tra_avg = tra_outlier';
  endif
end