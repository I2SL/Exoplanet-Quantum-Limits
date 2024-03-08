function fz_nm = FourierZernike(r,th,n,m)
    % Computes the Fourier Transform of the Zernike modes 
    % at the focal plane of a circular pupil
    %
    % r,th:     [N,1] vector of positions polar coordinates
    % n,m:      [1,M] vector of mode indices
    %
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024
    
    fz_nm = (-1).^(n/2 + abs(m)) .* sqrt(n+1) .* FZRadial(r,n) .* FZAzimuthal(th,m);
end

function z = FZRadial(r,n)
    % Computes the radial function of the Fourier Transformed Zernikes
    nn = repmat(n,[size(r,1),1,size(r,3)]);
    rr = repmat(r,[1,size(n,2),1]);

    % sinc-bessel in polar
    J = besselj(nn+1,2*pi*rr) ./ (sqrt(pi) * rr);
    
    % fill in singularities
    J(r==0 , n+1 == 1) = sqrt(pi);
    J(r==0 , n+1 > 1) = 0;
    
    % radial function
    z = J;
end

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

