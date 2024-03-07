load("CoronagraphEigenmodes.mat")
%{
%% display coronagraph eigenstates
nmodes = 5;
figure
t = tiledlayout(3,nmodes);
for k = 1:nmodes
    nexttile(k)
    phplot(PI_eig(:,:,k))
    title([sprintf('$\\chi_{%i}(\\mathbf{r})$',k-1)],'interpreter','latex')
    xticks([]); yticks([])
    axis square
    if k==1
        ylabel('PIAACMC','interpreter','latex')
    end
    
    nexttile(k+nmodes)
    phplot(VC_eig(:,:,k))
    xticks([]); yticks([])
    axis square
    if k==1
        ylabel('Vortex','interpreter','latex')
    end

    nexttile(k+2*nmodes)
    bar(abs([PI_lam(k),VC_lam(k);nan,nan]).^2)
    xlim([.5,1.5])
    xticks([1])
    xticklabels(sprintf('\\tau_{%i}',k-1))
    ylabel('Mode Transmission','interpreter','latex')
    if k==nmodes
        legend({'PIAACMC','Vortex'},'interpreter','latex')
    end
    %set(gca,'yscale','log')

end

plot([PI_lam.^2,VC_lam.^2]);

% plot the log of the mode intensities
figure 
colormap(slanCM('bubblegum'))
t = tiledlayout(3,nmodes);
for k = 1:nmodes
    nexttile(k)
    imagesc((abs(PI_eig(:,:,k)).^2))
    title([sprintf('$\\chi_{%i}(\\mathbf{r})$',k-1)],'interpreter','latex')
    xticks([]); yticks([])
    axis square
    if k==1
        ylabel('PIAACMC','interpreter','latex')
    end
    
    nexttile(k+nmodes)
    imagesc((abs(VC_eig(:,:,k)).^2))
    xticks([]); yticks([])
    axis square
    if k==1
        ylabel('Vortex','interpreter','latex')
    end

    nexttile(k+2*nmodes)
    bar(abs([PI_lam(k),VC_lam(k);nan,nan]).^2)
    xlim([.5,1.5])
    xticks([1])
    xticklabels(sprintf('\\tau_{%i}',k-1))
    ylabel('Mode Transmission','interpreter','latex')
    if k==nmodes
        legend({'PIAACMC','Vortex'},'interpreter','latex')
    end
    set(gca,'yscale','log')

end

%% plot eigenmode transmission coefficients
figure
t = tiledlayout(2,2);
nexttile(1)
%phplot(PI_U)
phplot(PI_C)
axis square

nexttile(2)
%phplot(VC_U)
phplot(VC_C)
axis square

nexttile(3)
PC_lam = ones(numel(PI_lam),1); PC_lam(n==0 & m==0) =0; 
bar(abs([PC_lam,PI_lam,VC_lam]).^2)
xticks(1:numel(PI_lam))
xticklabels(arrayfun(@(j) sprintf('\\chi_{%i}',j-1),1:numel(PI_lam),'UniformOutput',false))

xlabel('Coronagraph Eigenmode','interpreter','latex')
ylabel('Transmission Coefficient $|\tau_{k}|^2$','interpreter','latex')
leg = legend({'Perfect','PIAACMC','Vortex'},'interpreter','latex');
leg.Location = 'NorthWest';


% plot zernike transmission coefficients
nexttile(4)

tz_PC = ones(1,numel(lind)); tz_PC(np==0 & mp==0) =0; 
tz_PI = sum(abs(PI_C(:,lind)).^2,1); % zernike transmission probability through PIAACMC 
tz_VC = sum(abs(VC_C(:,lind)).^2,1); % zernike transmission probability through vortex

bar([tz_PC;tz_PI;tz_VC]')
ylim([0,1])
xticks(1:numel(lind))
xticklabels(arrayfun(@(j) sprintf('\\psi_{%i %i}',np(j),mp(j)),1:numel(lind),'UniformOutput',false))

xlabel('Zernike Mode','interpreter','latex')
ylabel('Zernike Mode Transmission','interpreter','latex')
leg = legend({'Perfect','PIAACMC','Vortex'},'interpreter','latex');
leg.Location = 'NorthWest';

%}

figure
nmodes = 4;
min_lam = 0;
id = 1:numel(PI_lam);
id_PI = abs(PI_lam).^2>=min_lam;
id_PI = id(id_PI);
id_VC = abs(VC_lam).^2>=min_lam;
id_VC = id(id_VC);


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

c = [
    [0.9290 0.6940 0.1250]
    [0.4660 0.8240 0.1880]
    [0 0.4470 0.8410]
    ];
nexttile(2*nmodes+1,[1,nmodes])
PC_lam = ones(numel(PI_lam),1); PC_lam(n==0 & m==0) =0; 
B = bar(abs([PC_lam,PI_lam,VC_lam]).^2);
%B(1).FaceColor = c(1,:);
%B(2).FaceColor = c(2,:);
%B(3).FaceColor = c(3,:);
xticks(1:numel(PI_lam))
xticklabels(arrayfun(@(j) sprintf('\\chi_{%i}',j-1),1:numel(PI_lam),'UniformOutput',false))
ylim([0,1.1])
xlabel('Coronagraph Eigenmode','interpreter','latex')
ylabel({'Eigenmode Transmission' ,'$|\tau_{k}|^2$'},'interpreter','latex')
leg = legend({'Perfect','PIAACMC','Vortex'},'interpreter','latex');
leg.Location = 'southeast';
%set(gca,'yscale','log')

%{
nexttile(3*nmodes+1,[1,nmodes])
np = [0,1,2,2,3,3];
mp = [0,1,0,2,1,3];
np = n;
mp = m;
lind = (np.*(np+2)+mp)/2 + 1;

tz_PC = ones(1,numel(lind)); tz_PC(np==0 & mp==0) =0; 
tz_PI = sum(abs(PI_C(:,lind)).^2,1); % zernike transmission probability through PIAACMC 
tz_VC = sum(abs(VC_C(:,lind)).^2,1); % zernike transmission probability through vortex


bar([tz_PC;tz_PI;tz_VC]')
ylim([0,1])
xticks(1:numel(lind))
xticklabels(arrayfun(@(j) sprintf('\\psi_{%i %i}',np(j),mp(j)),1:numel(lind),'UniformOutput',false))

xlabel('Zernike Mode','interpreter','latex')
ylabel('Zernike Mode Transmission','interpreter','latex')
leg = legend({'Perfect','PIAACMC','Vortex'},'interpreter','latex');
leg.Location = 'NorthWest';

%}



