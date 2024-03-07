% Visualize Fourier Zernike Modes

% rayleigh scaling parameter
rl = 1.22/2;

% Zernike Mode Indices
n_max = 4;
num_modes = (n_max+1)*(n_max+2)/2;
[n,m] = ZernikeIndices(n_max);

% coordinate grid
width = 10*rl;
n_grid = 501;
x = width*linspace(-.5,.5,n_grid);
[X,Y] = meshgrid(x);
[Th,R] = cart2pol(X,Y);

% FZ modes
FZ = FourierZernike(R(:),Th(:),n,m);
FZ = reshape(FZ,[n_grid,n_grid,num_modes]);

% plot
shading interp
colormap hot
dims = [n_max+1,2*n_max+1];
t = tiledlayout(dims(1),dims(2));
t.TileSpacing = 'compact';
t.Padding = 'compact';
for j = 1:num_modes
    k = sub2ind(fliplr(dims),m(j)+ n_max+1, n(j)+1);
    nexttile(k)
    %h = pcolor(abs(FZ(:,:,j)).^2);
    h = pcolor(log(abs(FZ(:,:,j)).^2));
    clim([-9,log(sqrt(pi))])
    title(sprintf('$(%i,%i)$',n(j),m(j)),'interpreter','latex')
    set(h, 'EdgeColor', 'none');
    grid off
    axis square
    axis off
end

shading interp
set(gcf,'renderer','painters');


function fz_nm = FourierZernike(r,th,n,m)
    fz_nm = (-1).^(n/2 + abs(m)) .* sqrt(n+1) .* FZRadial(r,n) .* FZAzimuthal(th,m);
end

function u = FZRadial(r,n)
    % Computes the radial function of the Fourier Transformed Zernikes
    nn = repmat(n,[size(r,1),1,size(r,3)]);
    rr = repmat(r,[1,size(n,2),1]);

    % sinc-bessel in polar
    J = besselj(nn+1,2*pi*rr) ./ (sqrt(pi) * rr);
    
    % fill in singularities
    J(r==0 , n+1 == 1) = sqrt(pi);
    J(r==0 , n+1 > 1) = 0;
    
    % radial function
    u = J;
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
