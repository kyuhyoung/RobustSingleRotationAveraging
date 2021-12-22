function transformation_averaging_from_text_file_of_x_y_z_rx_ry_rz(fn_tra_rot, mode_rot_avg)  
  %x = randn(1000,1);
  %edges = [-10 -2:0.25:2 10];
  %h = hist(x, edges,  "facecolor", "r", "edgecolor", "b" );
  %h = hist(x);
  %h
  [li_homo_mat, li_tra_meter, li_rot_vec] = load_translations_and_rotations_from_text_file_of_x_y_z_rx_ry_rz(fn_tra_rot);  
  disp('')
  rot_mat_avg = RobustSingleRotationAveraging_from_list_of_angle_axis_vectors(li_rot_vec, mode_rot_avg)
  rot_vec_avg = logarithm_map(rot_mat_avg)
  disp('')
  tra_avg_meter = median_translation_from_list_of_translation_vectors(li_tra_meter)
  li_mm_dif_wrt_avg_tra = 1000.0 * compute_list_of_distance_wrt_point(li_tra_meter, tra_avg_meter');
  li_deg_dif_wrt_avg_rot = compute_list_of_degree_difference_wrt_rotation(li_rot_vec, rot_vec_avg);
  edge_tra_mm = [0 : 1.0 : 8];
  cnt_mm_dif = histc(li_mm_dif_wrt_avg_tra, edge_tra_mm); prob_mm_dif = cnt_mm_dif / sum(cnt_mm_dif);
  edge_rot_deg = [0 : 0.2 : 1.6];
  cnt_deg_dif = histc(li_deg_dif_wrt_avg_rot, edge_rot_deg); prob_deg_dif = cnt_deg_dif / sum(cnt_deg_dif);
  %li_mm_dif_wrt_avg_tra'
  %prob_mm_dif'
  %li_deg_dif_wrt_avg_rot'
  %prob_deg_dif'
  %edge_tra_mm
  %prob_mm_dif
  figure
  bar(edge_tra_mm, prob_mm_dif, 'histc')
  xlabel('dist in mm')
  ylabel('prob')
  axis([-Inf Inf 0 1])
  set(gca, 'XTick', edge_tra_mm)
  grid on 
  
  figure
  bar(edge_rot_deg, prob_deg_dif, 'histc')
  xlabel('deg dif')
  ylabel('prob')
  axis([-Inf Inf 0 1])
  set(gca, 'XTick', edge_rot_deg)
  grid on 
  
endfunction


function li_deg_dif = compute_list_of_degree_difference_wrt_rotation(li_rot_vec_cell, rot_vec_wrt)
  n_rot = length(li_rot_vec_cell);
  li_deg_dif = -ones(n_rot, 1);
  rot_mat_wrt = RotationFromAngleAxis(rot_vec_wrt);
  for iR = 1 : n_rot
    rot_mat = RotationFromAngleAxis(li_rot_vec_cell{iR});
    rot_mat_dif = rot_mat * rot_mat_wrt';
    rot_vec_dif = logarithm_map(rot_mat_dif);
    li_deg_dif(iR) = norm(rot_vec_dif) * 180.0 / pi;  
  end
  %li_deg_dif
  %pause(100);
endfunction


function li_dist = compute_list_of_distance_wrt_point(li_p_cell, pnt)
  n_pt = length(li_p_cell);
  li_dist = -ones(n_pt, 1);
  for iP = 1 : n_pt
    li_dist(iP) = norm(li_p_cell{iP} - pnt);
  end   
  %li_dist
  %pause(100);
endfunction
