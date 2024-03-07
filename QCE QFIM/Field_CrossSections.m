
% Change interpreters to latex
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% constants
R = 1;
rl = 1.22/2;
u = R*linspace(-2,2,1001);
x = rl*linspace(-3,3,1001);

r_delta = rl;
b = 5e-2;

% pupil plane functions
phi_s = -2*pi*b*r_delta*u .* (abs(u)<=R);
phi_e = +2*pi*(1-b)*r_delta*u .* (abs(u)<=R);

% image plane functions
psi_s = besselj(1,2*pi*(x + b*r_delta)) ./ (sqrt(pi)*(x + b*r_delta));
psi_e = besselj(1,2*pi*(x - (1-b)*r_delta)) ./ (sqrt(pi)*(x - (1-b)*r_delta));


% plot 
tiledlayout(1,2)
nexttile(1)
plt = plot(u/R, phi_s,'-k', u/R, phi_e,':r');
hold on
yline(0)
hold off
xlim([-2,2])
plt(1).LineWidth = 2;
plt(2).LineWidth = 2;
yticks([])
xticks([-1,0,1])
xticklabels({'$-1$','$0$','$1$'})
axis square
xlabel({'$\mathbf{u}\cdot\bar{\mathbf{r}}_{\Delta}$','Cross Section Over Entrance Pupil'},'interpreter','latex')
ylabel({'Field Phase Profile'},'interpreter','latex')
leg = legend('Star','Planet');
leg.Location = 'northwest';

nexttile(2)
plt = plot(x/rl,psi_s,'-k',x/rl,psi_e,':r');
hold on
yline(0)
hold off
plt(1).LineWidth = 2;
plt(2).LineWidth = 2;
yticks([])
xticks([-2,-1,0,1,2])
xlim([-3,3])
axis square
xticklabels({'$-2\sigma$','$-\sigma$','$0$','$\sigma$','$2\sigma$'})
xlabel({'$\mathbf{r}\cdot\bar{\mathbf{r}}_{\Delta}$','Cross-Section Over Image Plane'},'interpreter','latex')
ylabel({'Field Amplitude Profile'},'interpreter','latex')
leg = legend('Star','Planet');
leg.Location = 'northwest';
leg.Location





