function R = RotationFromUnitAxisAngle(unit_axis, angle)
    
    if (angle==0)
        R = eye(3);
    else
        %unit_axis        
        so3 = SkewSymmetricMatrix(unit_axis);
        %so3
        %pause(100);        
        R = eye(3)+so3*sin(angle)+so3^2*(1-cos(angle));
    end
end



function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
end
