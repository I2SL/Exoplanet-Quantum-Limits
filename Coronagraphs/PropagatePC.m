function psi = PropagatePC(psi,X,Y)
    %----------------------------------------------------
    %-------------------- INPUTS ------------------------
    %----------------------------------------------------
    % psi    : input field(s) at focal plane [NxNxD]
    % X,Y    : meshgrid of focal plane coordinates [NxN]
    %----------------------------------------------------
    %-------------------- OUTPUTS -----------------------
    %----------------------------------------------------
    % psi : output field(s) on image plane


    % ----------------------------------------------------
    % dimensions
    [N,~,D] = size(psi);
    %[M,~] = size(Kx);

    % differentials
    d2x = (X(1,2) - X(1,1))^2;

    % cartesian to polar coordinates
    [Th,R]  = cart2pol(X,Y);
    % ----------------------------------------------------

    % get the PSF mode
    psf = FourierZernike(R(:),Th(:),0,0);

    % remove projection of the mode from the field(s)
    psi = reshape(psi,[N^2,D]);
    psi = psi - psf'*psi*d2x .* psf;

    % fourier transform to get pupil field
    %Psi = ctsIFT_2D(Kx(:), Ky(:), X(:), Y(:), d2x, psi);

    % reshape outputs
    psi = reshape(psi,[N,N,D]);
    %Psi = reshape(Psi,[M,M,D]);
end

