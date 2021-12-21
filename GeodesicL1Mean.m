function R = GeodesicL1Mean(R_input, b_outlier_rejection, n_iterations, thr_convergence)
    
    % 1. Initialize
    
    n_samples = length(R_input);
    
    vectors_total = zeros(9,n_samples);
    for i = 1:n_samples
        vectors_total(:,i)= R_input{i}(:);
    end
    s = median(vectors_total,2);
    %s = 
    %  0.873260
    %  0.014008
    %  -0.473210
    %  0.276570
    %  0.786950
    %  0.534230
    %  0.371340
    %  -0.605090
    %  0.692390    
    
    %vectors_total
    %s
    %pause(100);
    
    [U,~,V] = svd(reshape(s, [3 3]));
    R = U*V.';
    if (det(R) < 0)
        V(:,3) = -V(:,3);
        R = U*V.';
    end
    %R
    %pause(100);
    % 2. Optimize
    
    for j = 1:n_iterations

        vs = zeros(3,n_samples);
        v_norms = zeros(1,n_samples);
        for i = 1:n_samples
            v =  logarithm_map(R_input{i} * R');
            %R_tmp = R_input{i} * R';
            %v_tmp = logarithm_map(R_input{i} * R');          
            %R_tmp
            %v_tmp
            %v
            %pause(100);
            v_norm = norm(v);
            vs(:,i) = v;
            v_norms(i) = v_norm;
        end
        
        % Compute the inlier threshold (if we reject outliers).
        thr = inf;
        if (b_outlier_rejection)
            sorted_v_norms = sort(v_norms);
            v_norm_firstQ = sorted_v_norms(ceil(n_samples/4));
            %v_norms  
            %0.0307985   0.0057731   0.0401294   0.0716570   0.0369661   2.4650449   1.9336392
            %sorted_v_norms  
            %0.0057731   0.0307985   0.0369661   0.0401294   0.0716570   1.9336392   2.4650449
            %v_norm_firstQ
            #0.030798
            %pause(100);
            if (n_samples <= 50)
                thr = max(v_norm_firstQ, 1);

            else
                thr = max(v_norm_firstQ, 0.5);
            end
        end
        %thr
        %vs
        %v_norms
        %pause(100);
        step_num = 0;
        step_den = 0;

        for i = 1:n_samples
            v =  vs(:,i);
            v_norm = v_norms(i);
            if (v_norm > thr)
                continue;
            end
            %if (i > 1)
            %  step_num
            %  v
            %  v_norm
            %  v_over_v_norm = v / v_norm
            %  step_num_after = step_num + v/v_norm
            %  pause(100);       
            %endif
            step_num = step_num + v/v_norm;
            step_den = step_den + 1/v_norm;
        end
        %step_num
        %step_den
        %pause(100);
        delta = step_num/step_den;
        delta_angle = norm(delta);
        delta_axis = delta/delta_angle;
        
        R_delta = RotationFromUnitAxisAngle(delta_axis, delta_angle);
        R = R_delta * R;
        if (delta_angle < thr_convergence)
            %delta
            %delta_angle
            %R_delta
            %R
            break;
        end
    end
    %pause(100);   
end


function out = logarithm_map(in)
    cos_theta = (trace(in)-1)/2;
    sin_theta = sqrt(1-cos_theta^2);
    theta = acos(cos_theta);
    ln_R = theta/(2*sin_theta)*(in-in');
    out = [ln_R(3,2);ln_R(1,3);ln_R(2,1)];
end

