function [psi, Psi] = PropagateVC(psi,X,Y,Kx,Ky)
    %----------------------------------------------------
    %-------------------- INPUTS ------------------------
    %----------------------------------------------------
    % psi    : input field(s) at focal plane [NxNxD]
    % X,Y    : meshgrid of focal plane coordinates [NxN]
    % Kx,Ky  : meshgrid of pupil plane coordinates [MxM]
    %----------------------------------------------------
    %-------------------- OUTPUTS -----------------------
    %----------------------------------------------------
    % psi : output field(s) on image plane
    % Psi : output field(s) at secondary pupil (after Lyot stop)
    

    % ----------------------------------------------------
    % dimensions
    [N,~,D] = size(psi);
    [M,~] = size(Kx);
    d = size(psi);

    % differentials
    d2x = (X(1,2) - X(1,1))^2;
    d2k = (Kx(1,2) - Kx(1,1))^2;

    % cartesian to polar coordinates
    [T,R]  = cart2pol(X,Y);
    [~,Kr] = cart2pol(Kx,Ky);
    % ----------------------------------------------------

    % vortex phase mask
    charge = 2;
    vortex_phase = exp(1i*charge*T);

    % apply the vortex mask
    psi = psi .* vortex_phase;

    % fourier transform the field
    Psi = ctsIFT_2D(Kx(:), Ky(:), X(:), Y(:), d2x, reshape(psi,[N^2,D]));
    
    % apply Lyot stop
    lyot = double(Kr(:) <= 1);
    Psi = Psi .* lyot;

    % fourier transform the field back to the image plane
    psi = ctsIFT_2D(X(:),Y(:),Kx(:),Ky(:),d2k,Psi);
    
    
    %%%%%
    % Reference Field (removing residuals)
    psi0 = FourierZernike(R(:),T(:),0,0);
    psi0 = psi0 .* vortex_phase(:);
    Psi0 = ctsIFT_2D(Kx(:), Ky(:), X(:), Y(:), d2x, psi0);
    Psi0 = Psi0 .* lyot;
    psi0 = ctsIFT_2D(X(:), Y(:), Kx(:), Ky(:), d2k, Psi0);

    % Make unit vectors
    psi0 = psi0/(sqrt(d2x * (psi0'*psi0)));
    Psi0 = Psi0/(sqrt(d2k * (Psi0'*Psi0)));


    % remove residual field
    psi = psi - d2x * conj(psi0'*psi) .* psi0;
    Psi = Psi - d2k * conj(Psi0'*Psi) .* Psi0;
    %%%%%
    
    % reshape outputs
    psi = reshape(psi,[N,N,d(3:end)]);
    Psi = reshape(Psi,[M,M,d(3:end)]);

end