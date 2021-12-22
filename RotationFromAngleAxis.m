function R = RotationFromAngleAxis(angle_axis)
  R = RotationFromUnitAxisAngle(angle_axis / norm(angle_axis), norm(angle_axis));
end
