function z = Zernike(u,th,n,m)
    % Computes the zernike modes over a circular pupil of radius 1
    %
    % u,th:     [N,1] vector of positions polar coordinates
    % n,m:      [1,M] vector of mode indices
    %
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024

    % clip outside of support
    z = ZRadial(u,n,m) .* ZAzimuthal(th,m) .* (u <= 1);
end


function z = ZRadial(u,n,m)
% Computes the radial function of the Zernike modes
% 
% r:   [N,1] vector of radial positions
% n,m: [1,M] vector of indices

assert(numel(n)==numel(m))

z = zeros(numel(u),numel(n));

for i = 1:numel(n)
    for j = 0:((n(i)-abs(m(i)))/2)
        numerator = (-1)^j * sqrt(n(i)+1) * factorial(n(i)-j);
        denominator = factorial(j) * factorial((n(i)+m(i))/2 - j) * factorial((n(i)-m(i))/2 - j);
        z(:,i) = z(:,i) + numerator/denominator * u.^(n(i)-2*j);
    end    
end

z = z/sqrt(pi);


end


function v = ZAzimuthal(th,m)
    % Computes the angular function of the Zernike modes
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