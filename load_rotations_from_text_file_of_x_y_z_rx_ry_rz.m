function [li_homo_mat, li_tra, li_rot_vec] = load_translations_and_rotations_from_text_file_of_x_y_z_rx_ry_rz(fn_x_y_z_rx_ry_rz)
  li_homo_mat = cell();
  li_tra = cell();
  li_rot_vec = cell();
  fid = fopen(fn_x_y_z_rx_ry_rz);
  while 1
    tline = fgetl(fid)
    if ~ischar(tline)
      break      
    endif
    if("#" == tline(1))
      continue
    endif    
    %li_str = split(tline, [",", " "])    
    li_str = strsplit(tline, {",", " "})   
    %class(li_str) 
    if 6 ~= length(li_str)
      continue
    endif
    x_y_z_rx_ry_rx_cell = cellfun(@str2num, li_str, 'un', 0);
    x_y_z_rx_ry_rx_mat = cell2mat(x_y_z_rx_ry_rx_cell);
    li_tra{end + 1} = x_y_z_rx_ry_rx_mat(1:3)
    li_rot_vec{end + 1} = x_y_z_rx_ry_rx_mat(4:end)
    homo_mat = x_y_z_rx_ry_rz_2_homogeneous_transform(x_y_z_rx_ry_rx_mat)
    li_homo_mat{end + 1} = homo_mat;
    pause(100);
  endwhile
  fclose(fid);
  pause(100)
  li_li_str = readlines(fn_x_y_z_rx_ry_rz);
  li_li_str
  pause(100);  
endfunction
