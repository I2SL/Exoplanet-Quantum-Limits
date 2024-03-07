n_max = 42;
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