rl = 1.22/2;
r_delta = rl*linspace(1e-2,1e1,1000)';
%b = [.5,.4,.3,.2,10.^-(1:2:11)];
b = [.5,.4,.3,.2,.1,10^-2,10^-3];
%b = [10^-3,10^-5,10^-7];
delta_p = 1-2*b; 

% QFI radial component
QFI_rr = (1-delta_p.^2)*pi^2 - 4*(1-delta_p.^2).*delta_p.^2 ...
                            .*(besselj(2,2*pi*r_delta)./(r_delta)).^2; 

% QFI radial component high-contrast limit
QFI_rr_hc = 4*pi^2*b.*(1-4*(besselj(2,2*pi*r_delta)./(pi*r_delta)).^2);
QFI_rr_hc_normalized = 4*pi^2*(1-4*(besselj(2,2*pi*r_delta)./(pi*r_delta)).^2);


% QFI angular component
QFI_th = (1-delta_p.^2)*pi^2.*r_delta.^2;
QFI_th_normalized = r_delta.^2;

% Plot QFI Components
figure
%c1 = flipud(hot(numel(b)+10));

c1 = [    1.0000    0.5714         0
    1.0000    0.4286         0
    1.0000    0.2857         0
    1.0000    0.1429         0
    1.0000         0         0
    0.8571         0         0
    0.71           0         0];

% radial component
hold on
for k = 1:numel(b)
    %plot(r_delta/rl,QFI_rr(:,k)./((1-delta_p(k).^2)*pi^2),'LineWidth',1,'Color',c1(k+round(size(c1,1)*0.4286),:))
    plot(r_delta/rl,QFI_rr(:,k)./((1-delta_p(k).^2)*pi^2),'LineWidth',1,'Color',c1(k,:))
end
plot(r_delta/rl,QFI_rr_hc_normalized/(4*pi^2),'k','LineWidth',1)
hold off

legend_names = {'.5','.4','.3','.2','$10^{-1}$','$10^{-2}$','$<10^{-3}$'};

%legend_names = [legend_names,arrayfun(@(j) sprintf('$10^{%.2g}$',log10(b(j))),5:numel(b),'UniformOutput',false)];
legend_names = [legend_names,'High-Contrast Limit'];
leg = legend(legend_names,'interpreter','latex');
title(leg,'Relative Brightness $b$','interpreter','latex')
xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
ylabel({'Radial QFI','$\mathcal{K}^{(r\phi)}_{11}/(1-\kappa^2)\pi^2$'},'interpreter','latex')
title('Quantum Fisher Information','interpreter','latex')
set(gca,'XScale','log')
ylim([0,1.1])
grid on

%% 
% Plot QFI Heat Map
rl = 1.22/2;                                % rayleigh limit
r_delta = rl*10.^(linspace(-2,1,100));      % star-planet separation
b = 10.^-(linspace(3,11,100));              % relative brightness
[R,B] = meshgrid(r_delta,b);

% QFI radial component
Delta_P = 1-2*B;
QFI_rr = (1-Delta_P.^2)*pi^2 .* ( 1 - 4.*Delta_P.^2 ...
                            .*(besselj(2,2*pi*R)./(pi*R)).^2);

QFI_th = (1-Delta_P.^2)*pi^2 .* R.^2;

% Localization Uncertainty per photon Standard deviation
%N_ref = 1e14;
N_ref = 1;
Sigma_Loc = sqrt(1./N_ref.*(1./QFI_rr + R.^2./QFI_th));
Sigma_Loc_rel = 1./R .* Sigma_Loc;
Sigma_Loc_rel_normalized = Sigma_Loc_rel .* sqrt(N_ref*(1-Delta_P.^2)*pi^2);


min_rel_loc_err = .1;
Pho_Count = ((1/min_rel_loc_err) * Sigma_Loc_rel).^2;

figure
cmap = colormap(hot);

xImg = log10(linspace(min(r_delta)/rl, max(r_delta)/rl, size(QFI_rr, 2)));
yImg = log10(linspace(min(b), max(b), size(QFI_rr, 1)));
%image(xImg, yImg,flipud(log10(QFI_rr)), 'CDataMapping', 'scaled');
%image(xImg, yImg,flipud(log10(Sigma_Loc_rel)), 'CDataMapping', 'scaled');
%image(xImg, yImg,flipud(log10(Pho_Count)), 'CDataMapping', 'scaled');
image(xImg, yImg,flipud(log10(1./Pho_Count)), 'CDataMapping', 'scaled');
cbar = colorbar;
cbar_TickLabels = cellfun(@(label) sprintf(['$10^{',label(2:end),'}$']),cbar.TickLabels,'UniformOutput',false);
cbar.TickLabels = cbar_TickLabels;
cbar.TickLabelInterpreter = 'latex';
%ylabel(cbar,{'Relative Localization Error $\sigma_{loc}/r_{\Delta}$','Assuming $N=10^{14}$ Photons Collected'},'interpreter','latex')
ylabel(cbar,{'Photons Required','for $\sigma_{loc}/r_{\Delta} = 0.1$'},'interpreter','latex')
axis square
ax = gca;
ax.XTick = -2:1;
ax.XTickLabel = arrayfun(@(j)sprintf('$10^{%i}$',j),-2:1,'UniformOutput',false)';
ax.YTickLabel = arrayfun(@(j)sprintf('$10^{%i}$',j),-fliplr(3:11),'UniformOutput',false)';
ax.TickLabelInterpreter= 'latex';
ax.YDir = "normal";
xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
ylabel('Relative Brightness $b$','interpreter','latex')
title({'Exoplanet Localization'},'interpreter','latex')



% compute isocontours of the number of photons required to reach desired
% localization stdev
N = 10.^(5:2:13)';
%sigma_loc = rl/1e2;
%b_N = (4*pi^2*sigma_loc^2*N).^(-1).*(1 + 1./(1-4*(besselj(2,2*pi*r_delta)./(pi*r_delta)).^2));

rel_err = .1;
b_N = (4*pi^2*(rel_err * r_delta).^2.*N).^(-1).*(1 + 1./(1-4*(besselj(2,2*pi*r_delta)./(pi*r_delta)).^2));

c1 = gray(numel(N));
hold on
for k = 1:numel(N)
    plot(log10(r_delta/rl),log10(b_N(k,:)),'Color',c1(k,:),'LineWidth',2)
end
hold off

