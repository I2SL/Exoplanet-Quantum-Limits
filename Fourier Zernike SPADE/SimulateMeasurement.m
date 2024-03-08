function [mode_count,varargout] = SimulateMeasurement(mu_pho,p,varargin)
    % ----------------
    % Required Inputs:
    % ---------------
    % mu_pho            - mean photon number 
    % p                 - modal probability distribution (PMF)
    % 
    % ----------------
    % Optional Inputs:
    % ----------------
    % quant_eff         -  a scalar between 0 and 1 (inclusive) indicating the efficiency with which detector converts photons into electrons
    % poiss_flag        -   trigger for sampling the the number of collected photons from a poisson distribution
    % dark_lambda       -   poisson rate of dark current at each photon counter
    % read_noise_sigma  -   standard deviation for zero-mean gaussian random variable modelling the read noise in the photoelectric conversion process
    % crosstalk_mtx     -   symmetric matrix representing the leakage between modes. Columns and rows of matrix must sum to 1. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    %%%%%%%%%%%%%%%%% Parser %%%%%%%%%%%%%%%%%%%%%%%%%%

    % defaults constitute a noiseless detector
    default_quant_eff = 1;
    default_poiss_flag = 0;
    default_dark_lambda = 0;
    default_read_noise_sigma = 0;
    default_crosstalk_mtx = eye(numel(p));

    P = inputParser;
    
    addRequired(P,'mu_pho');
    addRequired(P,'p');
    
    addOptional(P,'quant_eff',default_quant_eff);
    addOptional(P,'poiss_flag',default_poiss_flag);
    addOptional(P,'dark_lambda',default_dark_lambda);
    addOptional(P,'read_noise_sigma',default_read_noise_sigma);
    addOptional(P,'crosstalk_mtx',default_crosstalk_mtx);
    
    parse(P, mu_pho, p, varargin{:});
    
    quant_eff = P.Results.quant_eff;
    poiss_flag = P.Results.poiss_flag;
    dark_lambda = P.Results.dark_lambda;
    read_noise_sigma = P.Results.read_noise_sigma;
    crosstalk_mtx = P.Results.crosstalk_mtx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % photon arrivals in each mode as poisson random variables with rates
    % mu_pho * mode_probability
    if poiss_flag
        mode_count = poissrnd(mu_pho*p);
    else
        mode_count = round(mu_pho*p);
    end

    % mix photon counts between modes according to crosstalk matrix
    assert(all(sum(crosstalk_mtx,2)-1 < 1e-9));
    mode_count = mode_count*crosstalk_mtx;

    % attenuate the number of registered photoelectrons by the quantum efficiency
    mode_count = round(quant_eff * mode_count);
    
    % add dark-current photoelectrons
    dc = poissrnd(dark_lambda, size(mode_count));
    mode_count = mode_count + dc;

    % add read noise photoelectrons
    rn = round(read_noise_sigma*randn(size(mode_count)));
    mode_count = mode_count + rn;
    mode_count(mode_count<0) = 0;
    
    % expand mode counts into measurement
    if nargout == 2
        varargout{1} = repelem(1:numel(p),mode_count);
    end
end