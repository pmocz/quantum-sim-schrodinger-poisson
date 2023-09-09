close all
clear all
clc

%% Make Figures for Paper
% Philip Mocz and Aaron Szasz (2020)


%% Carry out Classical Simulation
n = 2; Nt = 8;
n = 3; Nt = 20;
n = 4; Nt = 75;
n = 5; Nt = 300;
n = 6; Nt = 1200;
n = 8; Nt = 1200*4*4;

tFinal = 3; 

[psiC, VC] = classicalSim(n,tFinal,Nt);

N = 2^n;
xlin_ref = linspace(0,8,N+1);
xlin_ref = xlin_ref(1:end-1);
tlin_ref = linspace(0,3,Nt + 1);

[X, Y] = meshgrid(tlin_ref, xlin_ref);
V = abs(psiC).^2;


%%
fh = figure;
imagesc([0 8], [0 tFinal], abs(psiC).^2')
%imagesc([0 8], [0 tFinal], arg(psiC)')
set(gca,'ydir','normal')
%caxis([0 1])
colorbar


%%
%stop

%% Load Fake Quantum Simulation Data
rhoQ = cell(6,1);


for n = 2:6
    
    filename = ['output/qsim' num2str(n) '.hdf5'];
    
    rhoQ{n} = h5read(filename,'/rhoQ');
    
    figure;
    imagesc(rhoQ{n}')
    
end



%%  errors
errors = cell(6,1);
for n = 2:6
    N = 2^n;
    xlin = linspace(0,8,N+1);
    xlin = xlin(1:end-1);
    sz = size(rhoQ{n});
    tlin = linspace(0,tFinal,sz(2));
    
    [Xq, Yq] = meshgrid(tlin,xlin);
    dx = xlin(2) - xlin(1);
    dt = tlin(2) - tlin(1);
    dA = dx * dt;
    
    sol = interp2(X,Y,V,Xq,Yq);
    
    errors{n} = norm(rhoQ{n} - sol) ./ norm(sol);
    
end


%%
figure;
for n = 2:6
loglog(n,errors{n},'rx')
hold on
end
ns = 2:6;
plot(ns,4*ns.^-2)


%% stop

%% Make Final Figure
close all

fig = figure;

for n = 2:6
subplot(3,5,n-1+10)
imagesc([0 8], [0 tFinal], rhoQ{n}')
caxis([0 9])
set(gca,'ydir','normal')
axis square
%axis off
set(gca, 'XTick', [])
set(gca, 'YTick', [])
th = xlabel(['$n=' num2str(n) '$'],'interpreter','latex','fontsize',14);
end

subplot(3,5,2+[1:3 6:8])
imagesc([0 8], [0 tFinal], V')
caxis([0 9])
set(gca,'ydir','normal')
axis square
title('classical reference','interpreter','latex','fontsize',14)
xlabel('$x$','interpreter','latex','fontsize',12)
ylabel('$t$','interpreter','latex','fontsize',12)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '$|\psi|^2$';
set(colorTitleHandle ,'String',titleString,'interpreter','latex','fontsize',14);
set(gca, 'YTick', 0:3)

subplot(3,5,[1:2 6:7])
ns = 1:6;
loglog(ns,3.5*ns.^-2,'k--')
axis([1.5 6.5 0.1 0.8])
hold on
for n = 2:6
qq(n-1) = errors{n};
end
loglog(2:6,qq,'ro-','linewidth',2,'markersize',8)

pbaspect([1 1.7 1])
xlabel('$n$','interpreter','latex','fontsize',12)
ylabel('L2 error','interpreter','latex','fontsize',12)
text(2.4, 0.2,'$\propto n^{-2}$','interpreter','latex','fontsize',12);

set(gca, 'XTick', 2:6)
set(gca, 'YTick', [.1 .2 .4 .6 .8])

colormap(hot)


% save figure
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
%fig.PaperSize = [fig_pos(3) fig_pos(4)+.4];
%print(fig,'mockSim.eps','-depsc2')
