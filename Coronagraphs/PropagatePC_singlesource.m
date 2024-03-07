function psi = PropagatePC_singlesource(psi,X,Y,xs,ys)
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

    % cartesian to polar coordinates
    [Th,R]  = cart2pol(X,Y);

    % expansion coefficient for projection onto fundamental mode
    [ts,rs] = cart2pol(xs(:),ys(:));
    z0 = reshape(FourierZernike(rs,ts,0,0)/ sqrt(pi),size(xs));

    % remove component in fundamental mode
    psi0 = reshape(FourierZernike(R(:),Th(:),0,0),size(psi,[1,2]));
    psi = psi - z0.*psi0; 
end