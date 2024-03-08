function Psi = PIAACMCapodization(Psi)
% Applies the 0'th order Slepian prolate spheroidal function over a
% circular disk. The dimensionless parameter, c = pi D a / 2 
% where D is the aperture diameter and a is the mask diameter.
% It is found that complete on-axis extinction occurs for a=1.06/D
% such that c = pi*1.06 / 2.
%
% The implementation here follows from the relations given in appendix C of
% Soummer et. al. 2003, "Stellar Coronagraphy with prolate apodized
% circular apertures"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions of pupil field
ndim = size(Psi,1);

% constants
isodd = mod(ndim,2);
a = 1.06;
u = linspace(0,.5,round(ndim/2));
r = u.';
du = u(2)-u(1);
U = repmat(u,[numel(r),1]);

% Eigen-equation kernel
K0 = (a*r.*besselj(0,a*pi*u).*besselj(1,a*pi*r) - a*u.*besselj(0,a*pi*r) .* besselj(1,a*pi*u)) ./ (4*pi*(r.^2-u.^2));
K0(u==r) = a^2/8 * (besselj(0,a*pi*U(u==r)).^2 + besselj(1,a*pi*U(u==r)).^2);
M = (2*pi)^2 .* u .* K0 * du;
[U,lam] = eig(M,'vector');

% solve prolate spheroid eigen-equation with eigenvalue equal to 0.5
[~,idx] = min(abs(abs(lam)-.5));
pupil_apodization = U(:,idx);

%{
% UNCOMMENT TO PLOT THE RADIAL APODIZATION FUNCTION OVER THE PUPIL
plot([r; D/2; 3/4 * D ], [pupil_apodization/pupil_apodization(r==0);0;0],'k','LineWidth',1.5)
ylim([0,1.2])
xlim([0,.75])

xlabel('Radial Pupil Coordinate $u [D]$','interpreter','latex')
ylabel('PIAACMC Pupil Apodization')
%}


% 2D pupil apodization

% interpolate radial apodization function around 2pi radians
M = sqrt(r.^2+r'.^2); % dummy matrix of different radii from origin
M = reshape(M, [1, length(r)^2]); % make into 1D list
A = interp1(r, pupil_apodization, M,'spline',0); % interpolate onto 1D list
A = reshape(A, [length(r),length(r)]); % make back into array

% put into all four quadrants
A = [rot90(A(:,(1+isodd):end)); A];
A = [rot90(rot90(A(:,(1+isodd):end))) A];


% apply apodizer to pupil plane field while ensuring energy throughput
% preserved
energy = sum(conj(Psi).*Psi,[1,2]); % energy before apodization
Psi = A .* Psi;
Psi = Psi .* sqrt(energy ./ sum(conj(Psi).*Psi,[1,2])); % preserve energy after apodization

end