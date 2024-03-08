clear
% load example file
load('example.mat')

% load all mat files and put them in an object
num_MC = numel(min_sep_frac);% number of Monte-Carlo configurations

%{
for j = 1:num_MC
    fname = ['Exoplanet_MC_',num2str(signal_pho),'signal_pho_',num2str(min_sep_frac(j)),'min_sep_frac.mat'];
    load(fname)
    figure
    hold on
    scatter(exoplanet_xy_est(:,1),exoplanet_xy_est(:,2),'filled','b')
    xlim(rl*[-1.1,1.1])
    ylim(rl*[-1.1,1.1])

    dist = vecnorm(exoplanet_xy_est - xy_src(2,:),2,2);
    exoplanet_xy_est(dist>.1*rl,:) = [-exoplanet_xy_est(dist>.1*rl,1),exoplanet_xy_est(dist>.1*rl,2)];
    scatter(exoplanet_xy_est(:,1),exoplanet_xy_est(:,2),'filled','r')
    
    hold off
    close all
    clear dist
    save(fname)
end
%}

for j = 1:num_MC
    fname = ['Exoplanet_MC_',num2str(signal_pho),'signal_pho_',num2str(min_sep_frac(j)),'min_sep_frac.mat'];
    load(fname)
    MC(j).exoplanet_xy_est = exoplanet_xy_est;
    MC(j).xy_src = xy_src;
end


% preprocess the estimates
threshold_distance = rl*.15;
for j = 1:num_MC
    temp = MC(j).exoplanet_xy_est;
    dist = vecnorm(MC(j).exoplanet_xy_est - MC(j).xy_src(2,:),2,2);
    temp = temp(dist<threshold_distance,:);
    MC(j).exoplanet_xy_est = temp;
end



%% VISUALIZATION
fig = figure;
set(fig,'Position', [488.0000  261.8000  791.4000  496.2000])

% plot the spiral path with localization error as the color
ax1 = axes;
legend_names_ax1 = {};

% spiral coordinates
[~,r_spiral_frac] = cart2pol(xy_spiral_frac(:,1),xy_spiral_frac(:,2));
r_spiral = rl * r_spiral_frac;

% QFI terms
delta_p = 1-2*exoplanet_brightness;
QFI_rr = (1-delta_p.^2)*pi^2 .* ( 1 - 4.*delta_p.^2 ...
                            .*(besselj(2,2*pi*r_spiral)./(pi*r_spiral)).^2);
QFI_th = (1-delta_p.^2)*pi^2 .* r_spiral.^2;

% localization error
sigma_loc = (sqrt(1./mean_pho.*(1./QFI_rr + r_spiral.^2./QFI_th)));
sigma_loc_rl = sigma_loc/rl;


hold on
scatter(ax1,xy_spiral_frac(:,1),xy_spiral_frac(:,2),150*sigma_loc_rl,sigma_loc_rl); legend_names_ax1{2} = '';
xline(0); legend_names_ax1{3} = '';
yline(0); legend_names_ax1{4} = '';
hold off

% plot the estimates
view(2)
ax2 = axes;
legend_names_ax2 = {};
lnm=1;


% rayleigh spot limit
th = linspace(0,2*pi,1000)';
r = ones(size(th));
[x_rl,y_rl] = pol2cart(th,r);
hold on
plot(ax2,x_rl,y_rl,':k','LineWidth',1.5)
legend_names_ax2{lnm} = 'Diffraction Limit'; lnm=lnm+1;
scatter(ax2,0,0,100,'k','filled','pentagram')
legend_names_ax2{lnm} = 'Star'; lnm=lnm+1;
hold off


hold on
for j = 1:num_MC
    % get the unique estimates and their relative occurence counts
    [~,xy_unique,frac] = groupcounts(MC(j).exoplanet_xy_est);
    xy_unique = cell2mat(xy_unique);
    scatter(ax2,xy_unique(:,1)/rl,xy_unique(:,2)/rl,3,frac/max(frac),'filled');
    if j ==1
        legend_names_ax2{lnm} = 'Exoplanet Estimates (SPADE)';
    else
        legend_names_ax2{lnm} = '';
    end
    lnm = lnm +1;

end
hold off


% linke the axes
linkaxes([ax1,ax2])
ax2.Visible = 'off';

% colormap adjustments
colormap(ax1,flipud(slanCM('winter')))
colormap(ax2,slanCM('ember'))
ax2.ColorScale='log';

% add legends
%leg1 = legend(ax1,legend_names_ax1);
leg2 = legend(ax2,legend_names_ax2,'interpreter','latex');

%leg1.Position = [0.2216 0.8895 0.1770 0.0371];
leg2.Position = [0.2212 0.8283 0.2696 0.0972];


% add colorbar for sigma and line up the axes
cb1 = colorbar(ax1,'Position',[.185 .11 .0375 .815]);
cb2 = colorbar(ax2,'Position',[.73 .11 .0375 .815]);
%clim(ax1,[min(sigma_loc_rl)-mean([min(sigma_loc_rl),max(sigma_loc_rl)])/10,max(sigma_loc_rl)]);
cb1.Ticks = cb1.Ticks([1,round(end/2),end]);
ylabel(cb1,{'Localization Uncertainty $\sigma_{loc} / \sigma$'},'interpreter','latex')
ylabel(cb2,{'Monte Carlo Estimate Density'},'interpreter','latex')
set([ax1,ax2],'Position',[.07 .11 .815 .815]);

% make axes square and add limits, labels
axis(ax1,'square')
axis(ax2,'square')
xlim([-1.03,1.03])
ylim([-1.03,1.03])
yticks(ax1,xticks);
xlabel(ax1,'$x/\sigma$','interpreter','latex')
%ylabel(ax1,'$y/\sigma$','interpreter','latex')


%% Plot of monte-carlo variances as a function of 
figure

% 2D rotation matrix
R = @(a) [cos(a) , -sin(a);
          sin(a) ,  cos(a) ];

xy_mu = zeros([num_MC,2]);
var_r = zeros([num_MC,1]);
var_r_perp = zeros([num_MC,1]);
sigma_loc_exp = zeros([num_MC,1]);

for j = 1:num_MC
    xy = MC(j).exoplanet_xy_est/rl;
    
    % fit to a gaussian
    gm = fitgmdist(xy,1);
    
    % get the major and minor axis variances corresponding to the source
    % position estimate
    xy_mu(j,:) = gm.mu;
    [tt_mu,~] = cart2pol(xy_mu(j,1),xy_mu(j,2)); 

    Sigma_diag = R(tt_mu)*gm.Sigma*R(-tt_mu);
    var_r(j) = Sigma_diag(1,1);
    var_r_perp(j) = Sigma_diag(2,2);

    sigma_loc_exp(j) = sqrt(trace(gm.Sigma));
end
r_mu = sqrt(sum(xy_mu.^2,2));

hold on
colormap(flipud(slanCM('winter')))
scatter(r_spiral_frac,sigma_loc_rl,10,sigma_loc_rl)
plot(0,0,'blue','LineWidth',3)
scatter(r_mu,sigma_loc_exp,'filled','r')
hold off
xlim([0,1.1])
ylim([0.02,.03])
xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','LineWidth',1,'interpreter','latex')
ylabel('Localization Uncertainty $[\sigma]$','LineWidth',1,'interpreter','latex')
leg = legend({'','Quantum Limit $\sigma_{loc}/\sigma$','Monte-Carlo GaussFit $\sqrt{{Tr}(\Sigma)}/\sigma$'},'interpreter','latex');
%title(leg,'Localization Uncertainty','interpreter','latex')
axis square

%{
%cmap = flipud(colormap(winter));
cmap = flipud(colormap(ax2,winter));
sigma_loc_wn = sigma_loc - min(sigma_loc);
sigma_loc_wn = (size(cmap,1)-1)*sigma_loc_wn /max(sigma_loc_wn);
step = 20;
hold on
for i = 1:step:(size(xy_spiral_frac,1)-step)
    plot(ax2,xy_spiral_frac(i:i+step,1),xy_spiral_frac(i:i+step,2),'LineWidth',2,'Color',cmap(round(sigma_loc_wn(i))+1,:))
end
hold off
%}

%leg = legend(legend_names);
%title(leg,'$r_{\Delta}/\sigma$','interpreter','latex')




%{
hold on
patch([xy_spiral_frac(:,1);xy_spiral_frac(:,1)],[xy_spiral_frac(:,2);(xy_spiral_frac(:,2)-1e-2)],[sigma_loc;sigma_loc],'edgecolor','interp')
patch(xy_spiral_frac(:,1),xy_spiral_frac(:,2),sigma_loc,'edgecolor','interp');
colormap("jet")
legend_names{2*j+2} = '';
hold off
%}


%{
% Add precision map as backdrop
n_grid = 201;
x = rl*linspace(-1,1,n_grid);
[X,Y] = meshgrid(x);
[Th,R] = cart2pol(X,Y);

% scene setup
exoplanet_brightness = 1e-9;

%  localization uncertainty
var_r = 1./(4*exoplanet_brightness*pi^2 * (1 - 4*(besselj(2,2*pi*R)./(pi*R)).^2));
var_th = 1./(4*exoplanet_brightness*pi^2 * R.^2);
var_loc = var_r + R.^2.*var_th;

sigma_loc = sqrt(var_loc);
frac_sigma_loc = sigma_loc ./ R;

colormap gray
%imagesc(1./frac_sigma_loc);
h = image([-1,1],[1,-1],min(sigma_loc(:))./sigma_loc);
set(gca,'YDir','normal');
cbar = colorbar;
title(cbar,'$1/\sigma_{loc}$','interpreter','latex')
uistack(h,'bottom')
%}





%{

% FIGURE 1: Fit Gaussian Mixture Model Estimated Exoplanet Locations
GMModel = fitgmdist(exoplanet_xy_est/rl,2);
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
figure
y = exoplanet_xy_est(:,1)<0;
hold on
gscatter(exoplanet_xy_est(:,1)/rl,exoplanet_xy_est(:,2)/rl,y);
scatter(xy_src(1,1),xy_src(1,2),'filled','k');
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
legend('True Planet Estimates','Degenerate Planet Estimates','Star','Gaussian Fit')
hold off
title('{\bf Monte-Carlo Exoplanet Localization Estimates}')

xlim([-1,1])
ylim([-1,1])
axis square
xlabel('$x/\sigma$','interpreter','latex')
ylabel('$y/\sigma$','interpreter','latex')
xticks([-.5,-.25, 0, .25,.5])
yticks([-.5,-.25, 0, .25,.5])
xticklabels({'-1/2','-1/4','0','1/4','1/2'})
yticklabels({'-1/2','-1/4','0','1/4','1/2'})

% FIGURE 2: Histogram of estimates
figure
i_xy = round(exoplanet_xy_est/(X(1,2)-X(1,1)) + ceil(size(X,1)/2));
[C,ia,~] = unique(i_xy,'rows');
Z = zeros(size(X));
Z(sub2ind(size(Z),C(:,2),C(:,1))) = ia;
[M,c] = contourf(X,Y,Z);
%}