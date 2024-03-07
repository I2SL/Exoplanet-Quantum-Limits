clear

% setup the star-planet system parameters
rl = 1.22*pi/(2*pi); % Rayleigh length
dr = rl/1000;        % step size in star-planet separation   
r_delta = (rl*.01):dr:(rl*2); % star-planet separation
r_delta = permute(r_delta,[3,1,2]);
phi_delta = pi/4;   % star-planet angle
b = 1e-6;

% compute SLD mode probabilities
P = SLDModeProb(r_delta, phi_delta, b);


%{
logP  = log(P);

dlogP = diff(logP,1,3)./dr;

CFI_per_mode = (dlogP).^2 .*P(:,:,1:end-1);

CFI = sum(CFI_per_mode,2);
%}


% compute SLD mode CFI
dP = diff(P,1,3)./dr;
CFI_per_mode = (dP.^2) ./ (( P(1,:,1:end-1)+P(1,:,2:end) ) / 2);
CFI = sum(CFI_per_mode,2);


% QFI
delta_p = 1 - 2*b;
QFI_rr = (1-delta_p.^2)*pi^2 - 4*(1-delta_p.^2).*delta_p.^2 ...
      .*(besselj(2,2*pi*r_delta)./(r_delta)).^2;



% show probabilities by mode
figure
scatter(repmat(squeeze(r_delta),[1,size(P,2)])/rl,squeeze(P)')
plot(squeeze(r_delta)/rl,squeeze(P)')
xlabel('$r_{\Delta}/\sigma$','interpreter','latex')
ylabel('$P_i$','interpreter','latex')
title({'Probability of Photon Detection by SLD mode',...
    ['$b:$ ',sprintf('%.2e',b)]},'interpreter','latex')
leg = legend({'1','2','3','4'});
title(leg,'SLD mode')


% show probability derivative by mode
figure
scatter(repmat(squeeze(r_delta(1:end-1)),[1,size(P,2)])/rl,squeeze(dP)')
xlabel('$r_{\Delta}/\sigma$','interpreter','latex')
ylabel('$\partial_{r_{\Delta}} P_i$','interpreter','latex')
title({'Derivative of Photon Detection Probability by SLD mode',...
    ['$b:$ ',sprintf('%.2e',b)]},'interpreter','latex')
leg = legend({'1','2','3','4'});
title(leg,'SLD mode')

% show CFI
figure;
%plot(squeeze(r_delta(1:end-1))/rl,log(squeeze(CFI)./((1-delta_p.^2)*pi^2)),'red',...
%     squeeze(r_delta)/rl,log(squeeze(QFI_rr)./((1-delta_p.^2)*pi^2)),'blue')
plot(squeeze(r_delta(1:end-1))/rl,(squeeze(CFI)./((1-delta_p.^2)*pi^2)),'red',...
     squeeze(r_delta)/rl,(squeeze(QFI_rr)./((1-delta_p.^2)*pi^2)),'blue')
xlabel('$r_{\Delta}/\sigma$','interpreter','latex')
ylabel('Fisher Information','interpreter','latex')
title({'Fisher Information for SLD mode',...
    ['$b:$ ',sprintf('%.2e',b)]},'interpreter','latex')
legend({'SLD CFI','QFI'})


figure
scatter(repmat(squeeze(r_delta),[1,size(P,2)])/rl,squeeze(P)')
xlabel('$r_{\Delta}/\sigma$','interpreter','latex')
ylabel('$P_i$','interpreter','latex')
title({'Probability of Photon Detection by SLD mode',...
    ['$b:$ ',sprintf('%.2e',b)]},'interpreter','latex')
leg = legend({'1','2','3','4'});
title(leg,'SLD mode')


W = GramMatrix(r_delta);
SLD_mat = POVM(W,b);


v= VideoWriter('Figures/SLD_eigenmodes','MPEG-4');
v.FrameRate = 100;
open(v);

figure
for k = 1:numel(r_delta)
    % show the SLD modes
    Show_SLD_Modes(r_delta(k),phi_delta,b,rl,SLD_mat(:,:,k));
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v)



function P = SLDModeProb(r_delta, phi_delta, b)

    % star and planet locations
    xs = -b   .* r_delta.*cos(phi_delta);  
    ys = -b   .* r_delta.*sin(phi_delta);
    xe = (1-b).* r_delta.*cos(phi_delta);  
    ye = (1-b).* r_delta.*sin(phi_delta);
    x = [xs; xe]; y = [ys; ye];

    % evaluate SLD mode functions at star and planet locations
    SLD_xy = SLD(x,y,r_delta,phi_delta,b);

    % probability per mode for an off-axis source
    % is proporational to eignmodes squared evaluated at shift location.
    P_per_source = abs(SLD_xy).^2 / pi;

    % total probability per SLD eigenmode
    P = pagemtimes([(1-b),b], P_per_source);

    dr = r_delta(2) - r_delta(1);
    dP_per_source = diff(P_per_source,1,3)./dr;
    dP = (1-b)*dP_per_source(1,:,:) + b*dP_per_source(2,:,:);


    CFI_per_mode = (dP.^2) ./ (( P(1,:,1:end-1)+P(1,:,2:end) ) / 2);
end

function SLD_xy = SLD(X,Y,r_delta,phi_delta,b)
    
    % Evaluate span functions at the coordinates X,Y
    Xs = X + b*(r_delta.*cos(phi_delta)); 
    Ys = Y + b.*(r_delta.*sin(phi_delta));
    Xe = X - (1-b).*(r_delta.*cos(phi_delta)); 
    Ye = Y - (1-b).*(r_delta.*sin(phi_delta));

    [TH,R] = cart2pol(X,Y);
    [~,Rs] = cart2pol(Xs,Ys);
    [~,Re] = cart2pol(Xe,Ye);

    psi_s = besselj(1,2*pi*Rs)./(sqrt(pi)*Rs);     psi_s(Rs == 0) = sqrt(pi);
    psi_e = besselj(1,2*pi*Re)./(sqrt(pi)*Re);     psi_e(Re == 0) = sqrt(pi);
    
    dpsi_s = -2*besselj(2,2*pi*Rs)./( sqrt(pi)*Rs.^2);      dpsi_s(Rs == 0) = (pi)^(3/2);
    dpsi_s = dpsi_s .* (b*r_delta + R.*cos(TH-phi_delta));
    
    dpsi_e = -2*besselj(2,2*pi*Re)./( sqrt(pi)*Re.^2);      dpsi_e(Re == 0) = (pi)^(3/2);
    dpsi_e = dpsi_e .* ((1-b)*r_delta - R.*cos(TH-phi_delta));
    
    % position space non-ortho wavefunctions evaluated at X,Y
    rspace_wavefns = [psi_s,psi_e,dpsi_s,dpsi_e];
    
    % Gram Matrix
    W = GramMatrix(r_delta);

    % SLD wavefunctions evaluated at X,Y
    SLD_mat = POVM(W,b);
    SLD_xy = pagemtimes(rspace_wavefns, SLD_mat); 
end



function W = GramMatrix(r_delta)
    % Analytic Gram Matrix
    delta = besselj(1,2*pi*r_delta) ./ (pi*r_delta); 
    delta_p = - ( besselj(1,2*pi*r_delta) - 3 .* besselj(3,2*pi*r_delta) ) ./ (pi*r_delta);
    h = -2 * besselj(2,2*pi*r_delta) ./ (pi*r_delta);
    
    % fill in undefined values at zero separation
    delta(r_delta == 0) = 1;
    delta_p(r_delta == 0 ) = -1;
    h(r_delta == 0 ) = 0;



    W = zeros(4,4,numel(r_delta));
    W(1,2,:) = delta;
    W(3,4,:) = delta_p;
    W(1,4,:) = h;
    W(2,3,:) = h;
    W = W + eye(4) + pagectranspose(W);
end


function SLD_weights = POVM(W,b)
    % closed form eigenvectors of Gram-Matrix
    a = W(1,2,:);
    d = W(1,4,:);
    c = W(3,4,:);

    w = sqrt( a.^2 + 4*d.^2 - 2.*a.*c + c.^2 );
    uno = ones([1,1,size(W,3)]);
    
    % eigen decomposition of gram matrix
    Vw = pagetranspose([...
                 (-a + c - w)./(2*d),  (a - c + w)./(2*d), - uno, uno;
                 ( a - c - w)./(2*d) , (a - c - w)./(2*d),   uno, uno;
                 (-a + c + w)./(2*d) , (a - c - w)./(2*d), - uno, uno;
                 ( a - c + w)./(2*d) , (a - c + w)./(2*d),   uno, uno;
         ]);
    Vw = Vw ./ abs(vecnorm(Vw));

    % eigenvalues
    Dw  = [ -w - a - c + 2;
            -w + a + c + 2;
            +w - a - c + 2;
            +w + a + c + 2]./2;


    %% check

    % intermediate factors
    delta = a;
    gamma = sqrt(1-4*b.*(1-b).*(1-abs(delta).^2));
    eta1 = (1-2*b) + gamma;   
    eta2 = (1-2*b) - gamma;
    alpha1 = eta1 .* sqrt(1 ./ (eta1.^2 + (2*b*delta).^2 + (2.*delta).^2 .*b .* eta1 ));
    alpha2 = eta2 .* sqrt(1 ./ (eta2.^2 + (2*b*delta).^2 + (2.*delta).^2 .*b .* eta2 ));
    beta1 = 2*b.*delta.*sqrt(1 ./ (eta1.^2 + (2*b*delta).^2 + (2.*delta).^2 .*b .* eta1 ));
    beta2 = 2*b.*delta.*sqrt(1 ./ (eta2.^2 + (2*b*delta).^2 + (2.*delta).^2 .*b .* eta2 ));
    lambda1 = (1 + gamma)/2; 
    lambda2 = (1 - gamma)/2;
    
    % eigenvectors for image space (non-zero eigenvalues)
    v1 = [alpha1; beta1; 0*uno; 0*uno];
    v2 = [alpha2; beta2; 0*uno; 0*uno];
    

    check_ortho  = squeeze(pagemtimes(pagemtimes(pagectranspose(v1),W),v2));
    %%

    % density operator and its derivative
    rho = zeros(4); 
    rho(1,1) = (1-b); 
    rho(2,2) = b;
    
    
    d_rho = zeros(4); 
    d_rho(3,1) = 1; d_rho(1,3) = 1; 
    d_rho(2,4) = 1; d_rho(4,2) = 1;
    d_rho = pi*b*(1-b).* d_rho;


    % Phi matrix is representation of non-orthogonal state vectors in
    % eigebasis of Gram matrix
    Psi_phi = sqrt(Dw) .* pagectranspose(Vw);

    % represent density and its derivative in terms of the orthonormal
    % states from the diagonalized gram matrix
    rho_phi = pagemtimes(pagemtimes(Psi_phi, rho), pagectranspose(Psi_phi));
    d_rho_phi = pagemtimes(pagemtimes(Psi_phi, d_rho), pagectranspose(Psi_phi));
    
    % get orthonormalized eigenvectors of density
    [V_rho,D_rho] = eigenshuffle(rho_phi);
    D_rho = permute(D_rho,[1,3,2]);
    D_rho(3:4,:,:) = 0; % enforce eigenvalues beyond the span of the star and planet state span to be zero
    
    
    %% check 2
    v1_phi = pagemtimes(Psi_phi,v1);
    v2_phi = pagemtimes(Psi_phi,v2);
    
    figure
    tiledlayout(2,2)
    nexttile
    imagesc(squeeze(v1_phi));
    title('Analytic')
    colorbar
    nexttile
    imagesc(squeeze(v2_phi));
    title('Analytic')
    colorbar
    nexttile
    imagesc(squeeze(V_rho(:,1,:)));
    title('Numeric')
    colorbar
    nexttile
    imagesc(-squeeze(V_rho(:,2,:)));
    title('Numeric')
    colorbar


    V_rho(:,1:2,:) = [v1_phi,v2_phi];
    D_rho(1:2,1,:) = [lambda1;lambda2];
    %%





    % Gram schmidt additions to the null space of rho
    temp1 = rand(4,1);
    temp2 = rand(4,1);
    n1 = temp1 - pagemtimes(V_rho(:,1:2,:),pagectranspose(pagemtimes(pagectranspose(temp1),V_rho(:,1:2,:))));
    n1 = n1./vecnorm(n1,2,1);
    V_rho(:,3,:) = n1;
    
    n2 = temp2 - pagemtimes(V_rho(:,1:3,:),pagectranspose(pagemtimes(pagectranspose(temp2),V_rho(:,1:3,:))));
    n2 = n2./vecnorm(n2,2,1);
    V_rho(:,4,:) = n2;
    

    % get SLD matrix
    V_drho_V = pagemtimes(pagemtimes(pagectranspose(V_rho) , d_rho_phi) , V_rho);
    L_phi = zeros(size(rho_phi));
    
    D = zeros(size(L_phi));
    for i = 1:4
        for j =1:4

            D(i,j,:) = D_rho(i,1,:) + D_rho(j,1,:);

            L_phi = L_phi + 2 * (V_drho_V(i,j,:) ./ (D(i,j,:) + (D(i,j,:) == 0)) .* (D(i,j,:) ~= 0)) .* pagemtimes(V_rho(:,i,:), pagectranspose(V_rho(:,j,:)));

        end
    end

    % enforce SLD to be hermitian
    L_phi = (L_phi + pagectranspose(L_phi))/2; 
    
    % get orthonormalized SLD eigenvectors
    [V_sld, D_sld] = eigenshuffle(L_phi);
    D_sld = permute(D_sld,[1,3,2]);

    % show histogram of off-diagonals for V_sld'*V_sld to see if the
    % eigenstates are orthogonal
    II = pagemtimes(pagectranspose(V_sld),V_sld);
    off_diag = abs(II-eye(4));
    check_ortho = sum(off_diag,'all');
    figure
    histogram(off_diag);
    title('Total energy in off-diagonal elements of $V_{sld}^{\dagger}V_{sld}$','interpreter','latex')


    %{

    delta_p = 1-2*b;
    % Compute QFI from SLD matrix and density matrixs
    QFI = pagemtimes(rho_phi,pagemtimes(L_phi,L_phi));
    QFI = sum(QFI.*eye(4),[1,2]);

    % Compute CFI for SLD modes
    alt_CFI = sum(D_sld.^2 .* eye(4) .* pagemtimes(pagectranspose(V_sld), pagemtimes(rho_phi, V_sld)), [1,2]);
    
    figure
    t=tiledlayout(1,2);
    
    nexttile
    plot(r_delta/rl, squeeze(alt_CFI)/((1-delta_p.^2)*pi^2))
    title('SLD CFI calculated as $\sum_i \ell_i^{2} <\ell_i|\hat{\rho}|\ell_i> $','interpreter','latex')
    set(gca,'XScale','log')
    ylim([0,1])
    axis square
    
    nexttile
    plot(r_delta/rl,squeeze(QFI)/((1-delta_p.^2)*pi^2))
    title('QFI calculated as $Tr(\hat{\rho} \hat{L}^2)$','interpreter','latex')
    set(gca,'XScale','log')
    ylim([0,1])
    axis square

    title(t,['Relative Brightness: b=',sprintf('$10^{%i}$',log10(b))],'interpreter','latex')
    %}

    % return the weights of the SLD modes with respect to the
    % non-orthogonalized mode basis used in the Gram Matrix
    SLD_weights = pagemtimes(Vw, Dw.^(-.5) .* V_sld);

    
    % Compute QFI from SLD matrix and density matrixs
    QFI = pagemtimes(rho_phi,pagemtimes(L_phi,L_phi));
    QFI = sum(QFI.*eye(4),[1,2]); 
    delta_p = 1-2*b;
    figure;
    plot(squeeze(QFI)/((1-delta_p.^2)*pi^2))
    ylim([0,1])
    title('QFI calculated as $Tr(\hat{\rho} \hat{L}^2)$','interpreter','latex')
    set(gca,'XScale','log')

    % check SLD mode probabilities
    P_m = sum(conj(V_sld).* pagemtimes(rho_phi, V_sld),1);
    figure;
    plot(squeeze(P_m)')
    ylim([0,1])
end



function Show_SLD_Modes(r_delta, phi_delta, b, rl, SLD_matrix)
    
    % k-space starting wavefunctions
    n_grid = 501;
    [Ux,Uy] = meshgrid(linspace(-1,1,n_grid));
    [Th,U]  = cart2pol(Ux,Uy); 
    Pu = (U <= 1) / sqrt(pi);
    
    psi_s = Pu .* exp(-1i*2*pi*b*U.*r_delta.*cos(Th-phi_delta));
    psi_e = Pu .* exp(+1i*2*pi*(1-b)*U.*r_delta.*cos(Th-phi_delta));
    dpsi_s = (-1i*2*U.*cos(Th-phi_delta)) .* psi_s;
    dpsi_e = (+1i*2*U.*cos(Th-phi_delta)) .* psi_e;
    
    kspace_wavefns = [psi_s(:),psi_e(:),dpsi_s(:),dpsi_e(:)];
    
    
    % r-space starting wavefunctions
    x = 10*rl*linspace(-.5,.5,n_grid);
    [X,Y] = meshgrid(x);
    Xs = X + b*(r_delta*cos(phi_delta)); Ys = Y + b*(r_delta*sin(phi_delta));
    Xe = X - (1-b)*(r_delta*cos(phi_delta)); Ye = Y - (1-b)*(r_delta*sin(phi_delta));
    
    [TH,R] = cart2pol(X,Y);
    [~,Rs] = cart2pol(Xs,Ys);
    [~,Re] = cart2pol(Xe,Ye);
    
    psi_s = besselj(1,2*pi*Rs)./(sqrt(pi)*Rs);     psi_s(Rs == 0) = sqrt(pi);
    psi_e = besselj(1,2*pi*Re)./(sqrt(pi)*Re);     psi_e(Re == 0) = sqrt(pi);
    dpsi_s = -2*besselj(2,2*pi*Rs)./( sqrt(pi)*Rs.^2);      dpsi_s(Rs == 0) = (pi)^(3/2);
    dpsi_s = dpsi_s .* (b*r_delta + R.*cos(TH-phi_delta));
    dpsi_e = -2*besselj(2,2*pi*Re)./( sqrt(pi)*Re.^2);      dpsi_e(Re == 0) = (pi)^(3/2);
    dpsi_e = dpsi_e .* ((1-b)*r_delta - R.*cos(TH-phi_delta));
    
    rspace_wavefns = [psi_s(:),psi_e(:),dpsi_s(:),dpsi_e(:)];
%{
    t1 = tiledlayout(2,2);
    for i = 1:4
        nexttile(i)
        %imagesc(abs(rPOVM_modes(:,:,i)).^2 / pi)
        %colorbar;
        %clim([0,1]);
        phplot(reshape(rspace_wavefns(:,i),[n_grid,n_grid]));
        axis square
        xticks([])
        yticks([])
        xlabel('$x/\sigma$','interpreter','latex')
        ylabel('$y/\sigma$','interpreter','latex')
    end

%}
    % make SLD images
    W = GramMatrix(r_delta);

    % SLD eigenmodes as POVM 
    kPOVM_modes = kspace_wavefns * SLD_matrix;
    rPOVM_modes = rspace_wavefns * SLD_matrix;


    % ensure eigenmodes are orthonormal
    check_ortho_k = kPOVM_modes'*kPOVM_modes * (Ux(1,2)-Ux(1,1))^2;
    check_ortho_r = rPOVM_modes'*rPOVM_modes * (X(1,2)-X(1,1))^2;

    
    % reshape modes to be 2D
    kPOVM_modes = reshape(kPOVM_modes,[n_grid,n_grid,4]);
    rPOVM_modes = reshape(rPOVM_modes,[n_grid,n_grid,4]);
    
    t2 = tiledlayout(2,2);
    for i = 1:4
        nexttile(i)
        %imagesc(abs(rPOVM_modes(:,:,i)).^2 / pi)
        %colorbar;
        %clim([0,1]);
        phplot(rPOVM_modes(:,:,i))
        axis square
        xticks([])
        yticks([])
        xlabel('$x/\sigma$','interpreter','latex')
        ylabel('$y/\sigma$','interpreter','latex')

    end
    title(t2,{'SLD Modes at Image Plane',...
            ['$r_{\Delta}/\sigma:$ ',sprintf('%.2f',r_delta/rl)],...
        ['$b:$ ',sprintf('%.2e',b)]},'interpreter','latex')

    drawnow
end
