function psi = PropagatePC(psi,X,Y)
    % Propagates the fields psi through the perfect coronagraph.
    %----------------------------------------------------
    %-------------------- INPUTS ------------------------
    %----------------------------------------------------
    % psi    : input field(s) at focal plane [NxNxD]
    % X,Y    : meshgrid of focal plane coordinates [NxN]
    %----------------------------------------------------
    %-------------------- OUTPUTS -----------------------
    %----------------------------------------------------
    % psi : output field(s) on image plane
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ----------------------------------------------------
    % dimensions
    [N,~,D] = size(psi);

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

    % reshape outputs
    psi = reshape(psi,[N,N,D]);
end

