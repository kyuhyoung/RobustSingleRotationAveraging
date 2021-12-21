%clc;
close all; clear all;
pkg load statistics;

%% Function definitions

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

%t1 = trace(R_true*R_geodesic')
%t2 = t1 - 1
%t3 = t2 / 2
%t4 = acosd(t3)
%t6 = acos(t3)
%t5 = abs(t4)

error_GeodesicL1Mean = abs(acosd((trace(R_true*R_geodesic')-1)/2));
error_ChordalL1Mean = abs(acosd((trace(R_true*R_chordal')-1)/2));
%R_geodesic
%R_chordal
disp(['Error (geodesic L1 mean) = ', num2str(error_GeodesicL1Mean), ' deg, took ', num2str(time_geodesic*1000), ' ms'])
disp(['Error (chordal L1 mean) = ', num2str(error_ChordalL1Mean), ' deg, took ', num2str(time_chordal*1000), ' ms' ])


disp('')



