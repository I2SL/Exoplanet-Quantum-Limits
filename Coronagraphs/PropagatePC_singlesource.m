function psi = PropagatePC_singlesource(psi,X,Y,xs,ys)
    % Propagates the field induced by a single point source
    % through the perfect coronagraph.
    %----------------------------------------------------
    %-------------------- INPUTS ------------------------
    %----------------------------------------------------
    % psi    : input field(s) at focal plane [NxNxD]
    % X,Y    : meshgrid of focal plane coordinates [NxN]
    % xs,ys  : coordinates of point source [D]
    %----------------------------------------------------
    %-------------------- OUTPUTS -----------------------
    %----------------------------------------------------
    % psi : output field(s) on image plane
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % cartesian to polar coordinates
    [Th,R]  = cart2pol(X,Y);

    % expansion coefficient for projection onto fundamental mode
    [ts,rs] = cart2pol(xs(:),ys(:));
    z0 = reshape(FourierZernike(rs,ts,0,0)/ sqrt(pi),size(xs));

    % remove component in fundamental mode
    psi0 = reshape(FourierZernike(R(:),Th(:),0,0),size(psi,[1,2]));
    psi = psi - z0.*psi0; 
end