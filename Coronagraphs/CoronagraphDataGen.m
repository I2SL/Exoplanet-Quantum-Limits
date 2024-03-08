% Generate data for computing Classical Fisher Information Plots and
% Classical Chernoff Exponents of SPADE, Perfect Coronagraph, PIAACMC, and
% Vortex Coronagraph. It computes the direct imaging probability
% distributions over the image plane for different contrasts and
% star-planet separation vector.
% 
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024

% save filename
filename = 'CoronagraphData.mat';

% Image Space
rl = 1.22/2;                                        % rayliegh limit in image plane units
ndim = 401;                                         % image plane dimensionality
num_rl = 80;                                        % image plane extent in number of rayleigh widths 
[X,Y] = meshgrid(rl*num_rl*linspace(-.5,.5,ndim));  % image plane coordinates (cartesian)
[T,R] = cart2pol(X,Y);                              % image plane coordinates (polar) 
d2x = (X(1,2)-X(1,1))^2;                            % image plane differential area
E2 = d2x*sum(abs(FourierZernike(R(:),T(:),0,0)).^2); % total ensquared energy of PSF from discretization


% Pupil Space
R0 = 1;                                             % aperture radius
D = 2*R0;                                           % aperture diameter
[Kx,Ky] = meshgrid(linspace(-R0,R0,ndim));          % pupil plane coordinates (cartesian)
[Kt,Kr] = cart2pol(Kx,Ky);                          % pupil plane coordinates (polar)
d2k = (Kx(1,2)-Ky(1,1))^2;                          % pupil plane differential area

% star-planet relative brightness samples
b = [1e-3,1e-6,1e-9];                               % contrasts to sample
kappa = 1-2*b;

% star-planet separation vector samples (polar coordinates)
r_delta = rl*10.^linspace(-2,log10(3),150);         % radial separation samples
a_delta = [0,2*pi*1e-3];                            % angular separation samples  


% parameter sample space
[B,Rd,Ad] = ndgrid(b,r_delta,a_delta);
B = permute(B,[4,5,1,2,3]);
Rd = permute(Rd,[4,5,1,2,3]);
Ad = permute(Ad,[4,5,1,2,3]);

% exoplanet positions
xe = (1-B).*Rd.*cos(Ad);
ye = (1-B).*Rd.*sin(Ad);

% star positions
xs = -B.*Rd.*cos(Ad);
ys = -B.*Rd.*cos(Ad); 

% focal plane field induced by planet and star positions
[Te,Re] = cart2pol(gpuArray(X + xe),gpuArray(Y + ye));  
psi_e = reshape(FourierZernike(Re(:),Te(:),0,0),size(Re));
clear Te Re % clear some memory after constructing fields
[Ts,Rs] = cart2pol(gpuArray(X + xs),gpuArray(Y + ys)); 
psi_s = reshape(FourierZernike(Rs(:),Ts(:),0,0),size(Rs)); 
clear Ts Rs % clear some memory after constructing fields

% perfect coronagraph
psi_e_PC_out = PropagatePC_singlesource(psi_e,X,Y,xe,ye);
psi_s_PC_out = PropagatePC_singlesource(psi_s,X,Y,xs,ys);
psi_e_PC_out = gather(psi_e_PC_out);
psi_s_PC_out = gather(psi_s_PC_out);


% Vortex coronagraph
psi_e = gather(psi_e);
psi_s = gather(psi_s);
[psi_e_VC_out, ~] = PropagateVC(psi_e,X,Y,Kx,Ky);
[psi_s_VC_out, ~] = PropagateVC(psi_s,X,Y,Kx,Ky);

clear psi_e psi_s

% pupil plane field induced by planet and star positions
Psi_e =  exp(1i*2*pi * (Kx .* xe + Ky .* ye)) .* (Kr<=R0) / sqrt(pi*R0^2);
Psi_s =  exp(1i*2*pi * (Kx .* xs + Ky .* ys)) .* (Kr<=R0) / sqrt(pi*R0^2);

clear xe ye xs ys

% PIAACMC
[psi_e_PI_out, ~] = PropagatePIAACMC(Psi_e,X,Y,Kx,Ky);
[psi_s_PI_out, ~] = PropagatePIAACMC(Psi_s,X,Y,Kx,Ky);

clear Psi_e Psi_s

% direct imaging probability
P_PC = (1-B) .* abs(psi_s_PC_out).^2 + B .* abs(psi_e_PC_out).^2;
P_PI = (1-B) .* abs(psi_s_PI_out).^2 + B .* abs(psi_e_PI_out).^2;
P_VC = (1-B) .* abs(psi_s_VC_out).^2 + B .* abs(psi_e_VC_out).^2;

% save data
save(filename)

% Run scripts for plotting the CFI and CCE
CoronagraphCFIM
PlotCFI
PlotCCE
