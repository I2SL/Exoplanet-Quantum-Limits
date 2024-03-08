%% Exoplanet Localization with Spatial Mode Demultiplexing (SPADE)
% 
% This script simulates the localization of an exoplanet around a star 
% emitters in a sub-Rayleigh Field-of-View. The script assumes a 
% circular entrance pupil and measurements are photon counts in the Fourier
% Zernike spatial modes (truncated to a maximum order defined by the user). 
% Various real-world noise sources can also be modelled by changing the
% inputs to the SimulateMeasurement function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% save directory
save_dir = 'MLE_exoplanet_localizations';
mkdir(save_dir)

% Zernike SPADE Setup
n_max = 10;                              % max radial Fourier-Zernike mode order
[n,m] = ZernikeIndices(n_max);
num_modes = (n_max+1)*(n_max+2)/2;
prob_fn = @(x_src,y_src) FourierZernikeProbs(x_src,y_src,n,m);

% image space setup
rl = 1.22/2;                             % rayleigh width
n_grid = 201;
x = rl*linspace(-1.1,1.1,n_grid);
[X,Y] = meshgrid(x);

% scene setup
num_sources = 2;                                                       % number of sources in constellation (star + planet)
offset_angle = pi/4;                                        
exoplanet_brightness = 1e-9;    % relative exoplanet contrast
b_src = [1,exoplanet_brightness]'/(1+exoplanet_brightness);   % source brightnesses
[x_exo_frac,y_exo_frac] = UniformSpiralSampling(30,2,.1,offset_angle); % exoplanet cartesian coordinates
min_sep_frac = vecnorm(xy_exo_frac,2,2); % minimum star-planet separation for all sampled ground-truth exoplanet locations
xy_exo = rl*xy_exo_frac;                 % exoplanet coordinates

% dense sampling of spiral for plotting purposes
[x_spiral_frac,y_spiral_frac] = UniformSpiralSampling(10000,2,.1,offset_angle);
xy_spiral_frac = [x_spiral_frac,y_spiral_frac];

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


function mode_probs = FourierZernikeProbs(s_x,s_y,n,m)
    % computes the FZ mode photon detection probabilities induced by
    % incoherent sources at locations s_x, s_y for modes with index [n,m]
    [th,r] = cart2pol(s_x,s_y);
    mode_probs = abs(FourierZernike(r,th,n,m)).^2 / pi;
end