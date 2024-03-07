%% Exoplanet Localization with Spatial Mode Demultiplexing (SPADE)
% 
% This script simulates the localization of an exoplanet around a star 
% emitters in a sub-Rayleigh Field-of-View. The script assumes a 
% circular entrance pupil and measurements are photon counts in the Fourier
% Zernik spatial modes (truncated to a maximum order defined by the user). 
% Various real-world noise sources can also be modelled by changing the
% inputs to the SimulateMeasurement function.
%
% Author: Nico Deshler
% 11-29-2023

clear; 
close all;

% HG mode-sorter setup
rl = 1.22/2;                             % rayleigh limit (characteristic width of HG basis)
n_max = 10;                              % number of Fourier-Zernike modes
[n,m] = ZernikeIndices(n_max);
num_modes = (n_max+1)*(n_max+2)/2;
prob_fn = @(x_src,y_src) FourierZernikeProbs(x_src,y_src,n,m);
save_dir = 'filler_rotated_201_gridsamples_10_nmax_500MC_300pho';
mkdir(save_dir)

% image/object space setup
n_grid = 201;
x = rl*linspace(-1.1,1.1,n_grid);
[X,Y] = meshgrid(x);

% scene setup
exoplanet_brightness = 1e-9;
num_sources = 2;                % number of sources in constellation
centroid_aligned = 1;           % align constellation centroid to optical axis
offset_angle = pi/4;
[x_exo_frac,y_exo_frac] = UniformSpiralSampling(30,2,.1,offset_angle);


%%
%select subset near axis
select = [3,5,8,16,21,27];
x_exo_frac = x_exo_frac(select);
y_exo_frac = y_exo_frac(select);



% rotate 45 degrees 
xy_exo_frac = [x_exo_frac,y_exo_frac];
xy_exo_frac = xy_exo_frac* [cos(pi/4),sin(pi/4);
                            -sin(pi/4),cos(pi/4)];

%%
min_sep_frac = vecnorm(xy_exo_frac,2,2);
xy_exo = rl*xy_exo_frac;

[x_spiral_frac,y_spiral_frac] = UniformSpiralSampling(10000,2,.1,offset_angle);
xy_spiral_frac = [x_spiral_frac,y_spiral_frac];
b_src = [1,exoplanet_brightness]'/(1+exoplanet_brightness);   % source brightnesses


% expectation Maximization setup
n_em_max = 30;
brite_flag = 0;
exoplanet_flag = 1;

% measurement setup
signal_pho = 300;
mean_pho = signal_pho*(1/exoplanet_brightness);
read_noise_sigma = 0;
dark_lambda = 0;
poiss_flag = 1;
crosstalk_mtx = eye(num_modes,num_modes);
        
% Run monte-carlos for localization
trials = 500;


for j = 1:size(xy_exo,1)

% scene
xy_src = [0 , 0;...
         xy_exo(j,1),xy_exo(j,2)];
xy_src = xy_src - sum(b_src.*xy_src,2);  % align to centroid
scene = [xy_src,b_src];                  % scene columns are [x,y,b]

exoplanet_xy_est = zeros(trials,2);
err = zeros(trials,1);

for t=1:trials
    % Simulate measurement
    p_s = prob_fn(xy_src(:,1),xy_src(:,2));
    mode_probs = sum(b_src .* p_s,1);
    mode_counts = SimulateMeasurement(mean_pho,mode_probs,...
                                    'poiss_flag',poiss_flag,...
                                    'read_noise_sigma', read_noise_sigma,...
                                    'crosstalk_mtx',crosstalk_mtx,...
                                    'dark_lambda',dark_lambda);
    
    % run expectation maximization
    [~, x_est_trc, y_est_trc, loglike_trc, count] = EM(mode_counts,num_sources,xy_src,b_src,prob_fn,X,Y,rl,n_em_max,brite_flag,exoplanet_flag);
    xy_est = [x_est_trc(:,end),y_est_trc(:,end)];
    
    % compute the localization error
    err_trc = arrayfun(@(k) LocalizationError(xy_src, [x_est_trc(:,k),y_est_trc(:,k)]),1:count);
    
    % normalize error by minimum separation
    err_trc = err_trc / (rl*min_sep_frac(j));

    % collect results
    exoplanet_xy_est(t,:) = xy_est(2,:);
    err(t) = err_trc(end);

end
exoplanet_xy_est = exoplanet_xy_est(1:t,:);
err = err(1:t);

fname = fullfile(save_dir,['Exoplanet_MC_',num2str(signal_pho),'signal_pho_',num2str(min_sep_frac(j)),'min_sep_frac.mat']);
save(fname,'xy_src','exoplanet_xy_est','err','min_sep_frac','exoplanet_brightness','mean_pho','rl','signal_pho','xy_spiral_frac')

end



%% FUNCTIONS

function fz_nm = FourierZernike(r,th,n,m)
    fz_nm = (-1).^(n/2 + abs(m)) .* sqrt(n+1) .* FZRadial(r,n) .* FZAzimuthal(th,m);
end

function u = FZRadial(r,n)
    % Computes the radial function of the Fourier Transformed Zernikes
    nn = repmat(n,[size(r,1),1,size(r,3)]);
    rr = repmat(r,[1,size(n,2),1]);

    % sinc-bessel in polar
    J = besselj(nn+1,2*pi*rr) ./ (sqrt(pi) * rr);
    
    % fill in singularities
    J(r==0 , n+1 == 1) = sqrt(pi);
    J(r==0 , n+1 > 1) = 0;
    
    % radial function
    u = J;
end

function v = FZAzimuthal(th,m)
    % Computes the angular function of the Fourier Transformed Zernikes
    mm = repmat(m,[size(th,1),1]);
    tt = repmat(th,[1,size(m,2)]);

    v = zeros(size(tt));
    
    % angular function
    c = cos(abs(mm).*tt);
    s = sin(abs(mm).*tt);
    
    v(mm>0) = sqrt(2) * c(mm>0);
    v(mm==0)= 1;
    v(mm<0) = sqrt(2) * s(mm<0);

    
end

function mode_probs = FourierZernikeProbs(s_x,s_y,n,m)
    % computes the FZ mode photon detection probabilities induced by
    % incoherent sources at locations src_coords for modes with index [n,m]
    [th,r] = cart2pol(s_x,s_y);
    mode_probs = abs(FourierZernike(r,th,n,m)).^2 / pi;
end