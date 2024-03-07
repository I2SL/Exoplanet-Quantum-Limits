% Image Space
rl = 1.22/2;
ndim = 201;
num_rl = 12;
[X,Y] = meshgrid(rl*num_rl*linspace(-.5,.5,ndim));
[T,R] = cart2pol(X,Y); 
d2x = (X(1,2)-X(1,1))^2; 

% Pupil Space
R0 = 1;
D = 2*R0;
[Kx,Ky] = meshgrid(linspace(-R0,R0,ndim));
[Kt,Kr] = cart2pol(Kx,Ky);
d2k = (Kx(1,2)-Ky(1,1))^2;

% make the Zernike Modes in pupil and image space
nmax = 6;
[n,m] = ZernikeIndices(nmax);
Z = reshape(Zernike(Kr(:),Kt(:),n,m),[size(Kr),numel(n)]);
FZ = reshape(FourierZernike(R(:),T(:),n,m),[size(R),numel(n)]);

% send the zernike modes through the coronagraphs 
[~,PI_out] = PropagatePIAACMC(Z,X,Y,Kx,Ky);
[VC_out,~] = PropagateVC(FZ,X,Y,Kx,Ky);

% take inner product of output fields with input modes (Coronagraph
% operator in FZ mode basis)
PI_C = reshape(Z,[numel(Kr),numel(n)])'*reshape(PI_out,[numel(Kr),numel(n)])*d2k;
VC_C = reshape(FZ,[numel(R),numel(n)])'*reshape(VC_out,[numel(R),numel(n)])*d2x;

% threshold the coronagraph operators to account for numerical instability
% in mode orthogonality
PI_C = PI_C .* (abs(PI_C).^2 >= max(abs(PI_C(:))).^2/10);
VC_C = VC_C .* (abs(VC_C).^2 >= max(abs(VC_C(:))).^2/10);

% enforce hermitticity
PI_C = (PI_C + PI_C')/2; 
VC_C = (VC_C + VC_C')/2;

% diagonalize
[PI_U,PI_lam] = eig(PI_C,'vector');
[VC_U,VC_lam] = eig(VC_C,'vector');

% sort eigenvalues/vectors in ascending order
[~,id] = sort(abs(PI_lam).^2);
PI_lam = PI_lam(id);
PI_U = PI_U(:,id);

[~,id] = sort(abs(VC_lam).^2);
VC_lam = VC_lam(id);
VC_U = VC_U(:,id);

% get eigenfunctions
PI_eig = reshape(reshape(FZ,[numel(Kr),numel(n)])*PI_U,size(FZ));
VC_eig = reshape(reshape(FZ,[numel(Kr),numel(n)])*VC_U,size(FZ));


save('CoronagraphEigenmodes.mat','PI_eig','PI_lam','VC_eig','VC_lam','PI_C','VC_C','rl','num_rl','ndim','d2x','d2k','nmax')