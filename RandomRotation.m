function R = RandomRotation(max_angle_rad)
    %disp(['max_angle_rad : ', num2str(max_angle_rad)])
    %disp(['max_angle_rad : ', num2str(max_angle_rad)]);
    %disp('shit');
    unit_axis = rand(3,1)-0.5;
    unit_axis = unit_axis/norm(unit_axis);
    angle = rand * max_angle_rad;
    %unit_axis
    %angle
    R = RotationFromUnitAxisAngle(unit_axis, angle);
    %R
    %pause(100);
end