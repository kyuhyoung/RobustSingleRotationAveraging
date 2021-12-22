function homo_mat = x_y_z_rx_ry_rz_2_homogeneous_transform(x_y_z_rx_ry_rz)
  homo_mat = eye(4);
  x_y_z = x_y_z_rx_ry_rz(1:3);
  rx_ry_rz = x_y_z_rx_ry_rz(4:end);
  rad_rot = norm(rx_ry_rz);
  axis_rot = rx_ry_rz / rad_rot;
  rot_mat = RotationFromUnitAxisAngle(axis_rot, rad_rot);
  homo_mat(1:3, 1:3) = rot_mat;
  homo_mat(1:3, 4) = x_y_z;
  %x_y_z_rx_ry_rz
  %rot_mat
  %homo_mat
  %pause(100);         
endfunction
