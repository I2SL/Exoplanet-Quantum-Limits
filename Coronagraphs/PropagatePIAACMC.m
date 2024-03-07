function [psi, Psi] = PropagatePIAACMC(Psi,X,Y,Kx,Ky)
    %----------------------------------------------------
    %-------------------- INPUTS ------------------------
    %----------------------------------------------------
    % Psi    : input field(s) at pupil plane [MxMxD]
    % X,Y    : meshgrid of focal plane coordinates [NxN]
    % Kx,Ky  : meshgrid of pupil plane coordinates [MxM]
    %----------------------------------------------------
    %-------------------- OUTPUTS -----------------------
    %----------------------------------------------------
    % psi : output field(s) on image plane
    % Psi : output field(s) at secondary pupil (after Lyot stop)

    % ----------------------------------------------------
    % dimensions
    [N,~] = size(X);
    [M,~,D] = size(Psi);
    d = size(Psi);
    
    % differentials
    d2x = (X(1,2) - X(1,1))^2;
    d2k = (Kx(1,2) - Kx(1,1))^2;

    % cartesian to polar coordinates
    [~,R]  = cart2pol(X,Y);
    [Kt,Kr] = cart2pol(Kx,Ky);
    % ----------------------------------------------------

    % apply apodizer
    Psi = PIAACMCapodization(Psi);

    % propagate pupil field to mask plane
    psi = ctsIFT_2D(X(:),Y(:),Kx(:),Ky(:),d2k,reshape(Psi,M^2,D));
    
    % apply circular phase mask
    a = 1.06;
    mask = double(1 - 2*(R(:) <= a/4)); % Roddier-Roddier pi-phase
    psi = mask .* psi;
    
    % propagate to pupil
    Psi = ctsIFT_2D(Kx(:),Ky(:),X(:),Y(:),d2x,psi);

    % apply lyot stop
    lyot = double(Kr(:) <= 1);%.9761;
    Psi = Psi .* lyot;

    % propagate to image plane
    psi = ctsIFT_2D(X(:),Y(:),Kx(:),Ky(:),d2k,Psi);
    
    
    %%%%
    % Reference Field (removing residuals)
    Psi0 = PIAACMCapodization(reshape(Zernike(Kr(:),Kt(:),0,0),[M,M]));
    psi0 = ctsIFT_2D(X(:),Y(:),Kx(:),Ky(:),d2k,Psi0(:));  
    psi0 = mask .*  psi0;
    Psi0 = ctsIFT_2D(Kx(:),Ky(:),X(:),Y(:),d2x,psi0);
    Psi0 = Psi0 .* lyot;
    psi0 = ctsIFT_2D(X(:),Y(:),Kx(:),Ky(:),d2k,Psi0);
    
    % Make unit vectors
    psi0 = psi0/(sqrt(d2x * (psi0'*psi0)));
    Psi0 = Psi0/(sqrt(d2k * (Psi0'*Psi0)));

   
    % remove residual field
    psi = psi - d2x*psi0'*psi .* psi0;
    Psi = Psi - d2k*Psi0'*Psi .* Psi0;
    %%%%%
    
    
    % reshape outputs
    psi = reshape(psi,[N,N,d(3:end)]);
    Psi = reshape(Psi,[M,M,d(3:end)]);
end

