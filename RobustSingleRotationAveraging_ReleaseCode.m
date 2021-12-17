%clc;
close all; clear all;
pkg load statistics;

%% Function definitions

function out = logarithm_map(in)
    cos_theta = (trace(in)-1)/2;
    sin_theta = sqrt(1-cos_theta^2);
    theta = acos(cos_theta);
    ln_R = theta/(2*sin_theta)*(in-in');
    out = [ln_R(3,2);ln_R(1,3);ln_R(2,1)];
end


function out = SkewSymmetricMatrix(in)
    out=[0 -in(3) in(2) ; in(3) 0 -in(1) ; -in(2) in(1) 0 ];
end

addpath(pwd);
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

function R = ProjectOntoSO3(M)   
    [U,~,V] = svd(M);
    R = U*V.';
    if (det(R) < 0)
        V(:,3) = -V(:,3);
        R = U*V.';
    end
end


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
        R = R_delta*R;
        if (delta_angle < thr_convergence)
            break;
        end
    end
end


function R = ChordalL1Mean(R_input, b_outlier_rejection, n_iterations, thr_convergence)
    
    % 1. Initialize
    n_samples = length(R_input);
    
    vectors_total = zeros(9,n_samples);
    for i = 1:n_samples
        vectors_total(:,i)= R_input{i}(:);
    end      

    s = median(vectors_total,2);
                
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

        update_medvec = s-s_prev;
        if (norm(update_medvec) < thr_convergence)
            break;
        end

    end
    
    R = ProjectOntoSO3(reshape(s, [3 3]));
    %R
    %s_tmp = reshape(s, [3 3])
    %R_tmp = ProjectOntoSO3(s_tmp)
    %pause(100);
end


%%addpath(pwd);

rand('state', 0.00);
%% Example: Average 100 rotations (50 inliers, 50 outliers)

%n_inliers = 50; n_outliers = 50;
n_inliers = 5; n_outliers = 2;
inlier_noise_level = 5; %deg;
R_true = RandomRotation(pi); 
%R_true
% R_true = [  0.880924   0.277119   0.383637;
%              0.013324   0.795787  -0.605429;
%              -0.473070   0.538449   0.697336 ]
%
% 1. Create input rotaions:
n_samples = n_inliers + n_outliers;
R_samples = cell(1, n_samples);
                
for i = 1:n_samples
    if (i <= n_inliers)
        % Inliers: perturb by 5 deg.
        axis_perturb = rand(3,1)-0.5;
        axis_perturb = axis_perturb/norm(axis_perturb);
        %angle_perturb = normrnd(0,inlier_noise_level/180*pi); 
        angle_perturb = unifrnd(-inlier_noise_level / 180 * pi,  inlier_noise_level / 180 * pi); 
        R_perturb = RotationFromUnitAxisAngle(axis_perturb, angle_perturb);
        %axis_perturb
        % 0.037644
        % -0.317405
        % 0.947543        
        %angle_perturb
        % -0.034328
        %angle_axis_perturb = axis_perturb * angle_perturb  
        % -0.0012923
        % 0.0108960
        % -0.0325276      
        %R_perturb
        % 0.9994117   0.0325142   0.0109149
        % -0.0325283   0.9994702   0.0011148
        % -0.0108728  -0.0014692   0.9999398
        %pause(100);
        R_samples{i} = R_perturb * R_true;
    else
        % Outliers: completely random.
        R_samples{i} = RandomRotation(pi); 
    end   
end

R_samples_2 = cell(1, n_samples);
R_samples_2(1) = [  0.87568	0.30871	0.37134; 
                    -0.015865	0.78695	-0.61681;
                    -0.48264	0.53423	0.69401 ];                    
R_samples_2(2) = [  0.88084	0.27657	0.38424;
                    0.014008	0.79603	-0.60509;
                    -0.47321	0.53837	0.6973 ];
R_samples_2(3) = [  0.89546	0.2712	0.35299;
                    0.012714	0.77708	-0.62928;
                    -0.44496	0.56798	0.69239 ];
R_samples_2(4) = [  0.85608	0.27909	0.43502;
                    0.060279	0.78202	-0.62033;
                    -0.51332	0.55727	0.65265 ];
R_samples_2(5) = [  0.87326	0.26435	0.40932;
                    0.030012	0.80926	-0.58668;
                    -0.48634	0.52461	0.69876 ];
R_samples_2(6) = [  -0.33503	0.24321	-0.91028;
                    -0.061857	0.95835	0.27882;
                    0.94017	0.14972	-0.30603  ];
R_samples_2(7) = [  -0.073099	0.96268	0.26058;
                    0.9324	0.1587	-0.32472;
                    -0.35396	0.21923	-0.9092 ];

R_samples = R_samples_2;
% 2-a. Average them using Hartley's L1 geodesic method 
% (with our initialization and outlier rejection scheme):

b_outlier_rejection = true;
n_iterations = 10;
thr_convergence = 0.001;
tic;
R_geodesic = GeodesicL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
time_geodesic = toc;

% 2-b. Average them using our approximate L1 chordal method 
% (with our initialization and outlier rejection shceme)

b_outlier_rejection = true;
n_iterations = 10;
thr_convergence = 0.001;
tic;
R_chordal = ChordalL1Mean(R_samples, b_outlier_rejection, n_iterations, thr_convergence);
time_chordal = toc;


% 3. Evaluate the rotation error (deg):

error_GeodesicL1Mean = abs(acosd((trace(R_true*R_geodesic')-1)/2));
error_ChordalL1Mean = abs(acosd((trace(R_true*R_chordal')-1)/2));

disp(['Error (geodesic L1 mean) = ', num2str(error_GeodesicL1Mean), ' deg, took ', num2str(time_geodesic*1000), ' ms'])
disp(['Error (chordal L1 mean) = ', num2str(error_ChordalL1Mean), ' deg, took ', num2str(time_chordal*1000), ' ms' ])
disp('')



