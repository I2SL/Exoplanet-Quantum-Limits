gamma = @(r) besselj(1,2*pi*r) ./ (pi*r);   % autocorrelation function
rl = 1.22/2;                                % rayleigh limit
%b = 10.^-(2:1:11)';                         % relative brightness
%b = 10.^-[2,2.5,3,3.5,4,5:2:11]';
%b = 10.^-[2,2.5,3,3.5,4]';
b = 10.^-[3,3.5,4,5]';
r_delta = rl*linspace(1e-2,10,1000);      % star-planet separation
r_s = r_delta .* b ;                        % star radial distance from origin
r_e = r_delta .* (1-b);                     % planet radial distance from origin

% quantum chernoff exponent
xi_Q = - log((1-b) .* gamma(r_s).^2 + b .* gamma(r_e).^2);
xi_Q(xi_Q < 0) = 0; % correct numerical error

% approximation of chernoff exponent in high contrast regime
xi_Q_hc = b .* (1 - gamma(r_delta).^2);
xi_Q_hc_normalized = (1 - gamma(r_delta).^2); % no brightness dependence

%
set(gca,'XScale','log')
%set(gca,'YScale','log')


%  plot the Quantum chernoff exponent and the high-contrast approximation
figure
c1 = flipud(hot(numel(b)+10));
hold on
for j = 1:numel(b)
    plot(r_delta/rl,xi_Q(j,:)./b(j),'LineWidth',1,'Color',c1(j+round(size(c1,1)*0.4286),:))
end
plot(r_delta/rl,xi_Q_hc_normalized,'-k','LineWidth',1)
hold off
xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
ylabel({'Quantum Chernoff Exponent','$\xi_{Q}/b$'},'interpreter','latex')
title('Quantum Chernoff Exponent','interpreter','latex')
set(gca,'XScale','log')
%set(gca,'YScale','log')
grid on
legend_names = arrayfun(@(j) sprintf('$10^{%.2g}$',log10(b(j))),1:numel(b),'UniformOutput',false);
legend_names{end} = ['$<$ ',legend_names{end}];
legend_names = [legend_names,'High-Contrast Limit'];
leg = legend(legend_names,'interpreter','latex');
leg.Location =  'southeast';
title(leg,'Relative Brightness $b$','interpreter','latex')






% Make a contour plot showing number of photons required to reach 1/e 
rl = 1.22/2;                                % rayleigh limit
r_delta = rl*10.^(linspace(-2,1,100));      % star-planet separation
b = 10.^-(linspace(3,11,100));              % relative brightness
[R,B] = meshgrid(r_delta,b);
R_s = R .*  B ;                        % star radial distance from origin
R_e = R .* (1-B);                     % planet radial distance from origin
Xi_Q = - log((1-B) .* gamma(R_s).^2 + B .* gamma(R_e).^2);
Xi_Q = Xi_Q - min(Xi_Q(:)) + 1e-14;

Pho_Count = 7./Xi_Q;

N_ref = 1e5;
P_err_min = exp(-N_ref*Xi_Q);
P_err_min = P_err_min + 1e-14;


figure
%colormap(flipud(hot))
colormap(hot)
xImg = log10(linspace(min(r_delta)/rl, max(r_delta)/rl, size(Xi_Q, 2)));
yImg = log10(linspace(min(b), max(b), size(Xi_Q, 1)));
%image(xImg, yImg,flipud(log10(Xi_Q)), 'CDataMapping', 'scaled');
%image(xImg, yImg,flipud(log10(P_err_min)), 'CDataMapping', 'scaled');
%image(xImg, yImg,flipud(log10(Pho_Count)), 'CDataMapping', 'scaled');
image(xImg, yImg,flipud(log10(Xi_Q/7)), 'CDataMapping', 'scaled');
cbar = colorbar;
cbar_TickLabels = cellfun(@(label) sprintf(['$10^{',label(2:end),'}$']),cbar.TickLabels,'UniformOutput',false);
%cbar.TickLabels = fliplr(cbar_TickLabels);
cbar.TickLabels = cbar_TickLabels;
cbar.TickLabelInterpreter = 'latex';
%ylabel(cbar,'QCE $\xi_{Q}$','interpreter','latex')
ylabel(cbar,{'Photons Required', 'for $P_{e,min} \sim e^{-7} \approx 0.001$'},'interpreter','latex')
axis square
ax = gca;
ax.XTick = -2:1;
ax.XTickLabel = arrayfun(@(j)sprintf('$10^{%i}$',j),-2:1,'UniformOutput',false)';
ax.YTickLabel = arrayfun(@(j)sprintf('$10^{%i}$',j),-fliplr(3:11),'UniformOutput',false)';
ax.TickLabelInterpreter= 'latex';
ax.YDir = "normal";
xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
ylabel('Relative Brightness $b$','interpreter','latex')
%title({'Exoplanet Hypothesis Testing','Isocontours for 99.9\% Decision Confidence','$P_{e,min} \sim e^{-7} \approx .001$'},'interpreter','latex')
title({'Exoplanet Detection'},'interpreter','latex')




% compute isocontours of the number of photons required to reach desired detection confidence
confidence = .999;
p_err = 1-confidence;
exp_factor = abs(round(log(p_err)));

N = 10.^(5:2:13)';
b_N = (exp_factor./N) ./ (1 - gamma(r_delta).^2);
%c1 = winter(numel(N));
c1 = gray(numel(N));
hold on
for k = 1:numel(N)
    plot(log10(r_delta/rl),log10(b_N(k,:)),'Color',c1(k,:),'LineWidth',2)
    %plot(log10(r_delta/rl),log10(b_N(k,:)),'Color','r','LineWidth',2)
end
hold off

% Visual Magnitude to Photon Flux
lambda = 1290e-9;           % [m] 1290nm center wavelength
dlambda = lambda*.01;       % [m] 1% bandwidth filter
aperture_diameter = 6.5;    % [m] aperture diameter (next generation space telescope)
mx = 5.357;                 % Star Visual Magnitude (Proxima Centauri)

% zero magnitude star flux
J0 = 193.1;                                     % [photons/cm^2/s/angstrom]
J0 = J0 * dlambda * 1e10 * (1e2)^2;             % [photons/m^2/s]
pho_flux_m0 = J0 * pi*(aperture_diameter/2)^2;  % [photons/s]

% star flux
pho_flux = pho_flux_m0*10.^(mx/-2.5);           % [photons/s]
integration_time = N./pho_flux; % [s]

% alternate
%{
J0 = 1.589e-20;             % [erg / cm^2 / s / Hz]
J0 = J0 * 1e-7 * (1e2)^2;   % [W / m^2 / Hz]
J0 = J0 * 1e26;             % [Jansky] = [10^-26 W / m^2 / Hz]

pho_flux = VisualMagnitude_to_PhotonFlux(mx,'lambda',lambda,'dlambda',dlambda,'J0',J0); % [photons/sec/m^2]
integration_time = N./(pho_flux * pi * (aperture_diameter/2)^2 ); % [s]
 
%}


