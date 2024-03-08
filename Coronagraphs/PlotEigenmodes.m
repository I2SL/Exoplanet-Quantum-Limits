% Plots the eigenmodes of the PIAACMC and Vortex Coronagraphs
%
% Author(s): Nico Deshler, University of Arizona
% Affiliation(s): Wyant College of Optical Sciences, University of Arizona
% Date: March 7, 2024

% load data
load("CoronagraphEigenmodes.mat")


nmodes = 4;                         % number of eigenmodes to display
min_lam = 0;                        % threshold for minimum eigenvalues
id = 1:numel(PI_lam);
id_PI = abs(PI_lam).^2>=min_lam;
id_PI = id(id_PI);
id_VC = abs(VC_lam).^2>=min_lam;
id_VC = id(id_VC);

figure
t = tiledlayout(3,nmodes);
for k = 1:nmodes
    nexttile(k)
    if k~=nmodes
        phplot(PI_eig(:,:,id_PI(k)))
    else
       phplot(PI_eig(:,:,id_PI(k)),1,0,1)
    end
    title([sprintf('$\\chi_{%i}(\\mathbf{r})$',k-1)],'interpreter','latex')
    xticks([]); yticks([])
    axis square
    
    if k==1
        ylabel('PIAACMC','interpreter','latex')
    end
    
    nexttile(nmodes+k)
    if k~=nmodes
        phplot(VC_eig(:,:,id_VC(k)))
    else
       phplot(VC_eig(:,:,id_VC(k)),1,0,1)
    end
    xticks([]); yticks([])
    axis square
    if k==1
        ylabel('Vortex','interpreter','latex')
    end
end

nexttile(2*nmodes+1,[1,nmodes])
PC_lam = ones(numel(PI_lam),1); PC_lam(n==0 & m==0) =0; 
B = bar(abs([PC_lam,PI_lam,VC_lam]).^2);
xticks(1:numel(PI_lam))
xticklabels(arrayfun(@(j) sprintf('\\chi_{%i}',j-1),1:numel(PI_lam),'UniformOutput',false))
ylim([0,1.1])
xlabel('Coronagraph Eigenmode','interpreter','latex')
ylabel({'Eigenmode Transmission' ,'$|\tau_{k}|^2$'},'interpreter','latex')
leg = legend({'Perfect','PIAACMC','Vortex'},'interpreter','latex');
leg.Location = 'southeast';