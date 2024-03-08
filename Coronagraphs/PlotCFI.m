% Plots the Classical Fisher Information curves for all coronagraphs
%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024

% load the CFI data
load('CoronagraphCFI.mat')

% make figure
figure;
tiledlayout(2,numel(b))

% coronagraph line colors
c = [
    [0.7216 0.1216 0.2235]
    [0.9290 0.6940 0.1250]
    [0.4660 0.8240 0.1880]
    [0 0.4470 0.8410]
    ];

for i = 1:numel(b)
    

    % plot the normalized radial CFI curves
    nexttile(i)
    
    hold on
    plot(r_delta/rl, QFI_r(i,:)/(1-kappa(i)^2)/pi^2,'k','LineWidth',1.5)
    plot(r_delta(2:end)/rl, SP_CFIr(i,:)/(1-kappa(i)^2)/pi^2,'-.','Color',c(1,:),'LineWidth',2);    
    plot(r_delta(2:end)/rl, PC_CFIr(i,:)/(1-kappa(i)^2)/pi^2,'--','Color',c(2,:),'LineWidth',2);
    plot(r_delta(2:end)/rl, PI_CFIr(i,:)/(1-kappa(i)^2)/pi^2,'-','Color',c(3,:),'LineWidth',2);
    plot(r_delta(2:end)/rl, VC_CFIr(i,:)/(1-kappa(i)^2)/pi^2,'-','Color',c(4,:),'LineWidth',2);
    hold off
    
    set(gca,'XScale','log')
    grid on
    axis square
    xlim([min(r_delta/rl),max(r_delta/rl)])
    ylim([0,1.1])
    title(sprintf('$b =10^{%i}$',log10(b(i))),'interpreter','latex')
    xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
    if i == 1
        ylabel({'Radial FI', '$\mathcal{I}^{(r\phi)}_{11} /(1-\kappa^2)\pi^2$'},'interpreter','latex')
    end
    if i==numel(b)
        leg = legend({'QFI $\mathcal{K}^{(r\phi)}_{11}$','SPADE','Perfect','PIAACMC','Vortex'},'interpreter','latex');
        leg.Location = 'SouthWest';
    end



    % plot the normalized angular CFI curves
    nexttile(i+numel(b))
    hold on
    plot(r_delta/rl, QFI_a(i,:)/(1-kappa(i)^2)/pi^2,'k','LineWidth',1.5)
    plot(r_delta/rl, SP_CFIa(i,:)/(1-kappa(i)^2)/pi^2,'-.','Color',c(1,:),'LineWidth',2);
    plot(r_delta/rl, PC_CFIa(i,:)/(1-kappa(i)^2)/pi^2,'--','Color',c(2,:),'LineWidth',2);
    plot(r_delta/rl, PI_CFIa(i,:)/(1-kappa(i)^2)/pi^2,'-','Color',c(3,:),'LineWidth',2);
    plot(r_delta/rl, VC_CFIa(i,:)/(1-kappa(i)^2)/pi^2,'-','Color',c(4,:),'LineWidth',2);
    hold off

    set(gca,'XScale','log')
    set(gca,'YScale','log')
    grid on
    axis square
    xlim([min(r_delta/rl),max(r_delta/rl)])
    ylim([10^-10,10^2])
    yticks([10^-10,10^-5,10^0])
    title(sprintf('$b =10^{%i}$',log10(b(i))),'interpreter','latex')
    xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
    if i == 1
        ylabel({'Angular FI', '$\mathcal{I}^{(r\phi)}_{22} /(1-\kappa^2)\pi^2$'},'interpreter','latex')
    end
    if i==numel(b)
        leg = legend({'QFI $\mathcal{K}^{(r\phi)}_{22}$','SPADE','Perfect','PIAACMC','Vortex'},'interpreter','latex');
        leg.Location = 'SouthEast';
    end


end