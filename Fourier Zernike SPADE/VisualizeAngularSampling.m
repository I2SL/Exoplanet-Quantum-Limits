% Visualize angular sampling sensitivity for a truncated Zernike mode
% basis. For running a localization estimator on Zernike SPADE measurements
% a good rule of thumb is to pick truncation order given by:
%
% n_max >= ceil(pi/2/phi_min)
%
% where phi_min is the minimum absolute angle (in radians) that an
% exoplanet makes with either the x or y axes of the coordinate grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_max = 3;
omega = 2*pi*linspace(0,1,10000);
S = sin((1:n_max)' .* omega);
C = cos((1:n_max)' .* omega);
T = [S;C];

figure; 
plot(omega/2/pi,T.^2,'k')
hold on
plot(omega/2/pi,max(T.^2,[],1),'b','LineWidth',2)
hold off
xlabel('Angle $\phi/2\pi$','interpreter','latex')
ylabel('Angular Mode Probability','interpreter','latex')
title(['$max(n)=',num2str(n_max),'$ Zernike Mode Angular Sensitivity Map'],'interpreter','latex')