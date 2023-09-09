close all
clear all
clc
%% Philip Mocz and Aaron Szasz (2020)
% Towards Cosmological Simulations of Dark Matter on Quantum Computers

% Please cite:

% @article{mocz2021toward,
%   title={Toward cosmological simulations of dark matter on quantum computers},
%   author={Mocz, Philip and Szasz, Aaron},
%   journal={The Astrophysical Journal},
%   volume={910},
%   number={1},
%   pages={29},
%   year={2021},
%   publisher={IOP Publishing}
% }

% Carry out classical and proposed simulated-quantum computer simulation of
% the 1D Schrodinger-Poisson equation
% N = 2^n


n = 2; Nt = 8;
n = 3; Nt = 20;
n = 4; Nt = 75;
n = 5; Nt = 300;

tFinal = 3;


%% Carry out Classical Simulation
[psiC, VC] = classicalSim(n,tFinal,Nt);
psiC0 = psiC(:,1);
VC0 = VC(:,1);


%% Plot solution

fh = figure;
imagesc([0 8], [0 tFinal], abs(psiC).^2')
%imagesc([0 8], [0 tFinal], arg(psiC)')
set(gca,'ydir','normal')
%caxis([0 1])
colorbar


%% Carry out Fake Quantum Simulation
%close all

N = 2^n;
nlambda = 2*n*(n-1);
Lbox = 8;
dx = Lbox / N;
state0 = zeros(2^n,1);
state0(1) = 1;
ntrial = 200;  400;100;2;16;32; 
dt = tFinal / Nt;

% first, set up the IC
%lambda0 = 0*pi*ones(nlambda,1);
lambda0 = 0*pi*ones(nlambda/2,1);
psiQ = @(lambda) sqrt(N) * Upsi( lambda ) * state0;
psiQerror = @(lambda) norm( psiQ(lambda) - psiC0, 2 );
psiQerror_reduced = @(lambda) norm( psiQ([lambda; zeros(nlambda/2,1)]) - psiC0, 2 );

%mymin
rng(42);
lambdaQ = lambda0;
runningError = psiQerror_reduced(lambdaQ);
for trial = 1:ntrial
    trial
    [tmpQ, tmpError] = fminunc(psiQerror_reduced,lambda0);
    if tmpError < runningError
        lambdaQ = tmpQ;
        runningError = tmpError;
    end
    lambda0 = 2*pi*rand(nlambda/2,1);
end
lambdaQ = [lambdaQ; zeros(nlambda/2,1)];


%lambdaQ = fminunc(Qerror,lambda0);
%lambdaQ = fminsearch(Qerror,lambda0);
%lambdaQ = fminbnd(Qerror,0,2*pi);
% Qerror(lambdaQ)

figure;
plot(abs(psiC0).^2,'b','linewidth',1.5);
hold on
plot( abs(psiQ(lambdaQ)).^2,'r--','linewidth',3);
hold off

lambdaQ0 = lambdaQ;


%% Also solve for the inital potential
ntheta = n*(n-1) + 1;
theta0 = 0*pi*ones(ntheta,1);
thetaNorm = sqrt(sum(VC0.^2));
theta0(end) = thetaNorm;



VQ = @(theta) theta(end) * UV( theta ) * state0;
rhoQ = abs(psiQ(lambdaQ)).^2;
VQerror = @(theta) norm(  (circshift(VQ(theta),-1) + circshift(VQ(theta),1) - 2*VQ(theta)) / dx^2 - (rhoQ-1), 2 );

VQ_fixed = @(theta) thetaNorm * UV( theta ) * state0;
VQerror_fixed = @(theta) norm(  (circshift(VQ_fixed(theta),-1) + circshift(VQ_fixed(theta),1) - 2*VQ_fixed(theta)) / dx^2 - (rhoQ-1), 2 );

VQerror_forced = @(theta) norm( VQ_fixed(theta) - VC0, 2 );



%mymin
rng(42);
thetaQ = theta0;
runningError = VQerror_forced(thetaQ);
for trial = 1:ntrial
    trial
    tmpQ = fminunc(VQerror_forced,theta0);
    tmpQ(end) = thetaNorm;
    tmpError = VQerror_forced(tmpQ);
    if tmpError < runningError
        thetaQ = tmpQ;
        runningError = tmpError;
    end
    theta0 = 2*pi*rand(ntheta,1);
    theta0(end) = thetaNorm;
end

thetaQb = thetaQ; % save

thetaQ = fminunc(VQerror,thetaQ);
%thetaQ(end) = thetaNorm;

figure;
plot(VC0,'b','linewidth',1.5);
hold on
Vnorm = mean(VQ(thetaQ));
plot( VQ(thetaQ) - Vnorm,'r--','linewidth',3);
hold off

thetaQ0 = thetaQ;



%% make sure the ICs are good!

% stop


%% Now that we have a good initial guesses, start the simulation:
ntrial = 10;
%close all

% fourier space variables
fftw('planner','measure');
kx = (-N/2:N/2-1)' * (2*pi/Lbox);
kSq = fftshift(kx.^2);
clear kx

lambdaQ = mod(lambdaQ0, 2*pi);
thetaQ = thetaQ0;
%thetaQ(1:end-1) = mod(thetaQ0(1:end-1), 2*pi);


lb = zeros(size(lambdaQ));
ub = ones(size(lambdaQ))*2*pi;
ub(end) = 10 * thetaQ0(end);

lambdaQsav = zeros(nlambda,Nt+1);
thetaQsav = zeros(ntheta,Nt+1);
lambdaQsav(:,1) = lambdaQ;
thetaQsav(:,1) = thetaQ;

rng(42);
options = optimoptions('fminunc','Algorithm','quasi-newton','MaxIterations',800,'MaxFunctionEvaluations',800*length(lambdaQ),'OptimalityTolerance',1e-8);
optionsSA = optimoptions('simulannealbnd','HybridFcn',{@fminunc, options});


tic
for i = 1:Nt
    i
    
    % evolve
    prev = psiQ(lambdaQ);
    prevV = VQ(thetaQ);

    % kinetic - drift
    prev = fftn(prev);
    prev = exp(dt/2 * -1.i * kSq) .*prev;
    prev = ifftn(prev);
    
    %evolveQerror = @(lambda) norm( psiQ(lambda) - prev +  (1.i*dt)*prevV.*prev -  (1.i*dt/2)*(circshift(prev,-1) + circshift(prev,1) - 2*prev) / dx^2, 2 );
    evolveQerror = @(lambda) norm( psiQ(lambda) - prev +  (1.i*dt)*prevV.*prev, 2 );
    %evolveQerror = @(lambda) norm( psiQ(lambda) - exp(-1.i * dt * prevV).*prev, 2 );
    
    lambdaQ = fminunc(evolveQerror,lambdaQ,options);
    
    
    
    % update potential
    rhoQ = abs(psiQ(lambdaQ)).^2;
    VQerror = @(theta) norm(  (circshift(VQ(theta),-1) + circshift(VQ(theta),1) - 2*VQ(theta)) / dx^2 - (rhoQ-1), 2 );
    
    thetaQ = fminunc(VQerror,thetaQ,options);
    
    
    % save values
    
    lambdaQsav(:,i+1) = lambdaQ;
    thetaQsav(:,i+1) = thetaQ;
    
    %% plot
    figure(101);
    plot(abs(psiC(:,i+1)).^2,'b','linewidth',1.5);
    hold on
    plot( abs(psiQ(lambdaQ)).^2,'r','linewidth',1.5);
    %plot( abs(next).^2,'g','linewidth',1);
    plot( VC(:,i+1)  -1,'b','linewidth',1);
    plot( VQ(thetaQ) -1 - mean(VQ(thetaQ)),'r--','linewidth',1);
    %hold off
    
    pause(0.001)
    

    
end
toc



%% Prep Output

rhoQ = zeros(N,Nt+1);
for i = 1:Nt+1
    rhoQ(:,i) = abs(psiQ(lambdaQsav(:,i))).^2;
end

potQ = zeros(N,Nt+1);
for i = 1:Nt+1
    potQ(:,i) = VQ(thetaQsav(:,i));
end

potQ = potQ - mean(potQ);

%% Save Output
filename = ['output/qsim' num2str(n) '.hdf5'];

if exist(filename,'file')
    delete(filename);
end
h5create(filename,'/lambdaQsav',size(lambdaQsav));
h5create(filename,'/thetaQsav',size(thetaQsav));
h5create(filename,'/rhoQ',size(rhoQ));
h5create(filename,'/potQ',size(potQ));

h5write(filename,'/lambdaQsav',lambdaQsav);
h5write(filename,'/thetaQsav',thetaQsav);
h5write(filename,'/rhoQ',rhoQ);
h5write(filename,'/potQ',potQ);


%% Comparison Figure;


fh = figure;
subplot(1,2,1)
imagesc([0 8], [0 tFinal], abs(psiC).^2')
caxis([0 9])
set(gca,'ydir','normal')
subplot(1,2,2)
imagesc([0 8], [0 tFinal], rhoQ')
caxis([0 9])
set(gca,'ydir','normal')
colorbar


%%
fh = figure;
subplot(1,2,1)
imagesc([0 8], [0 tFinal], VC')
caxis([-4 4])
set(gca,'ydir','normal')
subplot(1,2,2)
imagesc([0 8], [0 tFinal], potQ')
caxis([-4 4])
set(gca,'ydir','normal')
colorbar
