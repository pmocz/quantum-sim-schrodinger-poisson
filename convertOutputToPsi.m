close all
clear all
clc

%% Save Psi
% Philip Mocz and Aaron Szasz (2020)


%% Carry out Classical Simulation
n = 2; Nt = 8;
n = 3; Nt = 20;
n = 4; Nt = 75;
n = 5; Nt = 300;
%n = 6; Nt = 1200;
%n = 8; Nt = 1200*4*4;

filename = ['output/qsim' num2str(n) '.hdf5'];
outfile = ['output/psiV' num2str(n) '.hdf5'];


tFinal = 3;

[psiC, VC] = classicalSim(n,tFinal,Nt);

N = 2^n;
nlambda = 2*n*(n-1);
Lbox = 8;
dx = Lbox / N;
state0 = zeros(2^n,1);
state0(1) = 1;
dt = tFinal / Nt;

% first, set up the IC
%lambda0 = 0*pi*ones(nlambda,1);
lambda0 = 0*pi*ones(nlambda/2,1);
psiQ = @(lambda) sqrt(N) * Upsi( lambda ) * state0;



%% Load Fake Quantum Simulation Data

V = h5read(filename,'/potQ');

lambdaQsav = h5read(filename,'/lambdaQsav');

psi = zeros(N,Nt+1);
for i = 1:Nt+1
    psi(:,i) = psiQ(lambdaQsav(:,i));
end



%% Quick Plot
figure;
subplot(2,1,1)
imagesc([0 8], [0 tFinal], [abs(psiC).^2'  abs(psi).^2'])
set(gca,'ydir','normal')
%caxis([0 1])
colorbar

subplot(2,1,2)
imagesc([0 8], [0 tFinal], [VC'  V'])
set(gca,'ydir','normal')
%caxis([0 1])
colorbar

%% save Results



h5create(outfile,'/psi',size(psi));
h5create(outfile,'/V',size(V));
h5create(outfile,'/psiC',size(psiC));
h5create(outfile,'/VC',size(VC));

h5write(outfile,'/psi',psi);
h5write(outfile,'/V',V);
h5write(outfile,'/psiC',psiC);
h5write(outfile,'/VC',VC);
