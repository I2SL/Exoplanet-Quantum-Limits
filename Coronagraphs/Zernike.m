function z = Zernike(r,th,n,m)
% compute the zernike modes over a circular support of radius 1
    
    % clip outside of support
    z = ZRadial(r,n,m) .* ZAzimuthal(th,m) .* (r <= 1);
end


function u = ZRadial(r,n,m)

assert(numel(n)==numel(m))

u = zeros(numel(r),numel(n));

for i = 1:numel(n)
    for j = 0:((n(i)-abs(m(i)))/2)
        numerator = (-1)^j * sqrt(n(i)+1) * factorial(n(i)-j);
        denominator = factorial(j) * factorial((n(i)+m(i))/2 - j) * factorial((n(i)-m(i))/2 - j);
        u(:,i) = u(:,i) + numerator/denominator * r.^(n(i)-2*j);
    end    
end

u = u/sqrt(pi);


end


function v = ZAzimuthal(th,m)
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