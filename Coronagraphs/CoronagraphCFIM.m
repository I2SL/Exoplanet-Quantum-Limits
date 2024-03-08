% Process coronagraph data to get radial and angular CFI. This script
% computes the numerical derivatives of the direct imaging probability
% distrubtions and calculates the CFI for the radial and angular estimation
% parameters as a function of their separation.
%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024

load('CoronagraphData.mat')


% get squared partial derivative w.r.t radius
PC_Dr2 = (diff(P_PC(:,:,:,:,1),1,4) ./ diff(Rd(:,:,:,:,1),1,4)).^2;
PI_Dr2 = (diff(P_PI(:,:,:,:,1),1,4) ./ diff(Rd(:,:,:,:,1),1,4)).^2;
VC_Dr2 = (diff(P_VC(:,:,:,:,1),1,4) ./ diff(Rd(:,:,:,:,1),1,4)).^2;

% get squared partial derivative w.r.t angle
PC_Da2 = (diff(P_PC,1,5) ./ diff(Ad,1,5)).^2;
PI_Da2 = (diff(P_PI,1,5) ./ diff(Ad,1,5)).^2;
VC_Da2 = (diff(P_VC,1,5) ./ diff(Ad,1,5)).^2;

% remove elements where PDF is zero
PC_Dr2 = PC_Dr2.*(P_PC(:,:,:,1:end-1,1) ~= 0); PC_Da2 = PC_Da2.*(P_PC(:,:,:,:,1:end-1) ~= 0);
PI_Dr2 = PI_Dr2.*(P_PI(:,:,:,1:end-1,1) ~= 0); PI_Da2 = PI_Da2.*(P_PI(:,:,:,:,1:end-1) ~= 0);
VC_Dr2 = VC_Dr2.*(P_VC(:,:,:,1:end-1,1) ~= 0); VC_Da2 = VC_Da2.*(P_VC(:,:,:,:,1:end-1) ~= 0);

%%%%%%% Radial CFI %%%%%%%%%%

%% PERFECT, PIAACMC, and VORTEX PROCESSING
% bias to first or second probability entry over which derivatives were calculated
lrp_w = 0.86; %0.935; 
lrp_w = permute([lrp_w;1-lrp_w],[5,4,3,2,1]);

P_PC_lerp = lrp_w.*P_PC(:,:,:,1:end-1,1) + (1-lrp_w).*P_PC(:,:,:,2:end,1);
P_PI_lerp = lrp_w.*P_PI(:,:,:,1:end-1,1) + (1-lrp_w).*P_PI(:,:,:,2:end,1);
P_VC_lerp = lrp_w.*P_VC(:,:,:,1:end-1,1) + (1-lrp_w).*P_VC(:,:,:,2:end,1);

PC_CFIrr = mean(PC_Dr2./P_PC_lerp,5);
PI_CFIrr = mean(PI_Dr2./P_PI_lerp,5);
VC_CFIrr = mean(VC_Dr2./P_VC_lerp,5);

% Angular CFI
PC_CFIaa = movmean(P_PC,2,5);
PC_CFIaa = PC_Da2./PC_CFIaa(:,:,:,:,2:end);
PI_CFIaa = movmean(P_PI,2,5);
PI_CFIaa = PI_Da2./PI_CFIaa(:,:,:,:,2:end);
VC_CFIaa = movmean(P_VC,2,5);
VC_CFIaa = VC_Da2./VC_CFIaa(:,:,:,:,2:end);

% remove nans induced by 0/0 error
PC_CFIrr(isnan(PC_CFIrr))=0; PC_CFIaa(isnan(PC_CFIaa))=0;
PI_CFIrr(isnan(PI_CFIrr))=0; PI_CFIaa(isnan(PI_CFIaa))=0;
VC_CFIrr(isnan(VC_CFIrr))=0; VC_CFIaa(isnan(VC_CFIaa))=0;


%% CFIs as function of separation

% get final CFIs                     
PC_CFIr = permute(d2x*sum(PC_CFIrr,[1,2]),[3,4,5,1,2]); PC_CFIa = permute(d2x*sum(PC_CFIaa,[1,2]),[3,5,4,1,2]);
PI_CFIr = permute(d2x*sum(PI_CFIrr,[1,2]),[3,4,5,1,2]); PI_CFIa = permute(d2x*sum(PI_CFIaa,[1,2]),[3,5,4,1,2]);
VC_CFIr = permute(d2x*sum(VC_CFIrr,[1,2]),[3,4,5,1,2]); VC_CFIa = permute(d2x*sum(VC_CFIaa,[1,2]),[3,5,4,1,2]);

PC_CFIr = medfilt1(PC_CFIr')';
PI_CFIr = medfilt1(PI_CFIr')';
VC_CFIr = medfilt1(VC_CFIr')';

PC_CFIa = squeeze(PC_CFIa);
PI_CFIa = squeeze(PI_CFIa);
VC_CFIa = squeeze(VC_CFIa);

%% SPADE CFI
nmax = 500;
[n,m] = ZernikeIndices(nmax);
th_delta = zeros(size(r_delta))+pi/4; 

% helper functions for computing Zernike CFIs
bsj = @(r,n) besselj(repmat(n,[size(r,1),1]), repmat(r,[1,size(n,2)]));
p_nm = @(r,th,n,m) abs(FourierZernike(r,th,n,m)).^2 / pi;
dp_nm_rr = @(r,th,n,m) 2 * bsj(2*pi*r,n+1) ./ repmat(r,[1,size(n,2)]) .* ( bsj(2*pi*r,n-1) - bsj(2*pi*r,n+3) ) .* (FZAzimuthal(th,m)).^2;
dp_nm_th = @(r,th,n,m) 4 * (bsj(2*pi*r,n+1) ./ (pi * repmat(r,[1,size(n,2)]))).^2 .* (n+1) .* abs(m) .* cos(abs(m).*th) .* sin(abs(m).*th) .* sign(-m);

% Compute CFI for each brightness ratio
for k = 1:numel(b)

    % star and planet coordinates
    r_s = b(k) * r_delta.';
    th_s = rem(th_delta.' + pi,2*pi);
    r_e = (1-b(k)) * r_delta.';
    th_e = th_delta.';

    % radial derivative of mode probabilities at star and planet locations
    dp_rr_s = dp_nm_rr(r_s,th_s,n,m);
    dp_rr_e = dp_nm_rr(r_e,th_e,n,m);

    % total radial derivative
    dP_nm_rr = b(k).*(1-b(k)).*(dp_rr_s + dp_rr_e);

    % azimuthal derivative of mode probabilities at star and planet
    % locations
    dp_th_s = dp_nm_th(r_s,th_s,n,m);
    dp_th_e = dp_nm_th(r_e,th_e,n,m);
    
    % total azimuthal derivative
    dP_nm_th = (1-b(k))*dp_th_s + b(k)*dp_th_e;
    
    % P_star
    Ps(:,:,k) = (1-b(k)).*p_nm(r_s,th_s,n,m);

    % P_exoplanet
    Pe(:,:,k) = b(k).*p_nm(r_e,th_e,n,m);

    % total probability
    P_nm(:,:,k) = Ps(:,:,k) + Pe(:,:,k);
    
    % CFI radial component
    CFI_nm_rr(:,:,k) = (dP_nm_rr).^2 ./ (P_nm(:,:,k) + realmin);
    CFI_rr(:,k) = sum(CFI_nm_rr(:,:,k),2);

    % CFI azimuthal component
    CFI_nm_th(:,:,k) = (dP_nm_th).^2 ./ (P_nm(:,:,k) + realmin);
    CFI_th(:,k) = sum(CFI_nm_th(:,:,k),2);
end

SP_CFIr = CFI_rr(2:end,:).';
SP_CFIa = CFI_th.';


%% QFIM
QFI_r = (1-kappa.'.^2)*pi^2 .* ((1) - (kappa.'.^2.*(2 * besselj(2, 2 *pi * r_delta) ./ (pi * r_delta) ).^2));
QFI_a = (1-kappa.'.^2)*pi^2 .* ((r_delta.^2));



%% SAVE
save('CoronagraphCFI.mat','QFI_r','QFI_a',...
                                 'SP_CFIr','SP_CFIa',...
                                 'PC_CFIr','PC_CFIa',...
                                 'PI_CFIr','PI_CFIa',...
                                 'VC_CFIr','VC_CFIa',...
                                 'r_delta','rl', ...
                                 'b','kappa','-v7.3')


%% FUNCTIONS
function v = FZAzimuthal(th,m)
    % Computes the angular function of the Fourier Transformed Zernikes
    mm = repmat(m,[size(th,1),1]);
    tt = repmat(th,[1,size(m,2)]);

    v = zeros(size(tt));
    
    % angular function
    c = cos(abs(mm).*tt);
    s = sin(abs(mm).*tt);
    
    v(mm>0) = sqrt(2) * c(mm>0);
    v(mm==0)= 1;
    v(mm<0) = sqrt(2) * s(mm<0);

end




