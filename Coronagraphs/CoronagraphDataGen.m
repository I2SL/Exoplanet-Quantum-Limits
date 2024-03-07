% Generates data for computing coronagraph CFIs and CCEs

% Image Space
rl = 1.22/2;
ndim = 401;
num_rl = 80;
[X,Y] = meshgrid(rl*num_rl*linspace(-.5,.5,ndim));
[T,R] = cart2pol(X,Y); 
d2x = (X(1,2)-X(1,1))^2; 
E2 = d2x*sum(abs(FourierZernike(R(:),T(:),0,0)).^2); % total ensquared energy from discretization


% Pupil Space
R0 = 1;
D = 2*R0;
[Kx,Ky] = meshgrid(linspace(-R0,R0,ndim));
[Kt,Kr] = cart2pol(Kx,Ky);
d2k = (Kx(1,2)-Ky(1,1))^2;

% star-planet relative brightness samples
b = [1e-3,1e-6,1e-9];
kappa = 1-2*b;

% star-planet separation vector samples (polar coordinates)
r_delta = rl*10.^linspace(-2,log10(3),150); % radial samples
a_delta = [0,2*pi*1e-3];                    % angular samples  


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

save(sprintf('CoronagraphData_FINAL2_%s.mat',datetime('now','Format',"yyyy-MM-dd-HH-mm-ss")))

CoronagraphCFIM
PlotCFI
