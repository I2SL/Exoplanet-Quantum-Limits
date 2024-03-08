% Visualize Fourier Zernike Modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Zernike Mode Indices
n_max = 4;                          % max radial order to plot
num_modes = (n_max+1)*(n_max+2)/2;
[n,m] = ZernikeIndices(n_max);

% Image Space
rl = 1.22/2;                                        % rayliegh limit in image plane units
ndim = 401;                                         % image plane dimensionality
num_rl = 6;                                        % image plane extent in number of rayleigh widths 
[X,Y] = meshgrid(rl*num_rl*linspace(-.5,.5,ndim));  % image plane coordinates (cartesian)
[Th,R] = cart2pol(X,Y);                              % image plane coordinates (polar) 

% FZ modes
FZ = FourierZernike(R(:),Th(:),n,m);
FZ = reshape(FZ,[ndim,ndim,num_modes]);

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
    h = pcolor(abs(FZ(:,:,j)).^2);
    title(sprintf('$(%i,%i)$',n(j),m(j)),'interpreter','latex')
    set(h, 'EdgeColor', 'none');
    grid off
    axis square
    axis off
end

shading interp
set(gcf,'renderer','painters');
