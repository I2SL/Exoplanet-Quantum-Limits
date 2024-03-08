% Plots the Classical Chernoff Exponent curves for all coronagraphs
%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024


% load data
load('CoronagraphData.mat','P_PI','P_VC','b','d2x','r_delta','rl')

% quantum chernoff exponent
gamma = @(r) besselj(1,2*pi*r) ./ (pi*r);   % autocorrelation function
r_e = (1-b').*r_delta;
r_s = b'.*r_delta;
xi_Q = - log((1-b') .* gamma(r_s).^2 + b' .* gamma(r_e).^2);
xi_Q(xi_Q < 0) = 0; % correct numerical error

% direct-imaging coronagraphs Classical Chernoff Exoponents
xi_SP = xi_Q;   % SPADE
xi_PC = xi_Q;   % Perfect Coronagraph
xi_PI = permute(-log(1-d2x*sum(P_PI(:,:,:,:,1),[1,2])),[3,4,1,2]);  % PIAACMC
xi_VC = permute(-log(1-d2x*sum(P_VC(:,:,:,:,1),[1,2])),[3,4,1,2]);  % Vortex


% coronagraph line colors
c = [
    [0.7216 0.1216 0.2235]
    [0.9290 0.6940 0.1250]
    [0.4660 0.8240 0.1880]
    [0 0.4470 0.8410]
    ];


figure
tiledlayout(1,numel(b))
for i = 1:numel(b)
    nexttile(i)
    hold on
    plot(r_delta/rl, xi_Q(i,:)/b(i),'k','LineWidth',2)
    plot(r_delta/rl, xi_SP(i,:)/b(i),'-.','Color',c(1,:),'LineWidth',2)
    plot(r_delta/rl, xi_PC(i,:)/b(i),'--','Color',c(2,:),'LineWidth',2)
    plot(r_delta/rl, xi_PI(i,:)/b(i),'-', 'Color',c(3,:),'LineWidth',2);
    plot(r_delta/rl, xi_VC(i,:)/b(i),'-', 'Color',c(4,:),'LineWidth',2);
    hold off

    title(sprintf('$b =10^{%i}$',log10(b(i))),'interpreter','latex')
    xlabel('Star-Planet Separation $r_{\Delta}/\sigma$','interpreter','latex')
    
    if i == 1
        ylabel({'Chernoff Exoponent','$\xi_{C}/b$'},'interpreter','latex')
    end

    if i==numel(b)
    leg = legend({'QCE $\xi_{Q}$','SPADE','Perfect','PIAACMC','Vortex'},'interpreter','latex');
    leg.Location = 'NorthWest';
    end

    
    ylim([0,1.1])
    xlim([min(r_delta/rl),max(r_delta/rl)])
    set(gca,'YScale','log')
    set(gca,'XScale','log')

    axis square
    grid on
end

