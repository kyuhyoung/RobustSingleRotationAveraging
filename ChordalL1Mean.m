function R = ChordalL1Mean(R_input, b_outlier_rejection, n_iterations, thr_convergence)
    
    % 1. Initialize
    n_samples = length(R_input);
    
    vectors_total = zeros(9,n_samples);
    for i = 1:n_samples
        vectors_total(:,i)= R_input{i}(:);
    end      

    s = median(vectors_total,2);
    %s
    %pause(100);
                
    % 2. Optimize
    for j = 1:n_iterations
        %vectors_total
        %s
        %vectors_total_minus_s = vectors_total - s
        %abs_vectors_total_minus_s = abs(vectors_total_minus_s)
        %sum_abs_vectors_total_minus_s = sum(abs_vectors_total_minus_s)
        %sum_abs_vectors_total_minus_s_0 = sum_abs_vectors_total_minus_s == 0
        %sum_sum_abs_vectors_total_minus_s_0 = sum(sum_abs_vectors_total_minus_s_0)
        %sum_sum_abs_vectors_total_minus_s_0_2 = sum(sum(abs(vectors_total - s))==0)
        if (sum(sum(abs(vectors_total - s))==0) ~= 0)
            %s
            s = s + rand(size(s, 1), 1) * 0.001;
            %s
            %pause(100);
        end

        v_norms = zeros(1,n_samples);
        for i = 1:n_samples
            v =  vectors_total(:,i)-s;
            v_norm = norm(v);
            v_norms(i) = v_norm;
        end

        % Compute the inlier threshold (if we reject outliers).
        thr = inf;
        if (b_outlier_rejection)
            sorted_v_norms = sort(v_norms);
            v_norm_firstQ = sorted_v_norms(ceil(n_samples/4));

            if (n_samples <= 50)
                thr = max(v_norm_firstQ, 1.356);
                % 2*sqrt(2)*sin(1/2) is approximately 1.356
            else
                thr = max(v_norm_firstQ, 0.7);
                % 2*sqrt(2)*sin(0.5/2) is approximately 0.7
            end
        end
        
        %thr
        %pause(100);

        step_num = 0;
        step_den = 0;

        for i = 1:n_samples
            v_norm = v_norms(i);
            if (v_norm > thr)
                continue;
            end
            step_num = step_num + vectors_total(:,i)/v_norm;
            step_den = step_den + 1/v_norm;
        end


        s_prev = s;
        s = step_num/step_den;
        %j
        %s
        %pause(100);
        
        update_medvec = s-s_prev;
        if (norm(update_medvec) < thr_convergence)
            break;
        end

    end
    %pause(100);
    
    R = ProjectOntoSO3(reshape(s, [3 3]));
    %R
    %s_tmp = reshape(s, [3 3])
    %R_tmp = ProjectOntoSO3(s_tmp)
    %pause(100);
end



function R = ProjectOntoSO3(M)   
    [U,~,V] = svd(M);
    R = U*V.';
    if (det(R) < 0)
        V(:,3) = -V(:,3);
        R = U*V.';
    end
end

