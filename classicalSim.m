function [psiOut, VOut] = classicalSim(n, tFinal, Nt)
%CLASSICALSIM simulate the 1D Schrodinger-Poisson system
%   using a spectral method
% Philip Mocz and Aaron Szasz (2020)

%% setup

N = 2^n;
Lbox = 8;


% spacing & grid
dx = Lbox / N;
x = (0:N-1)' * dx;   % Note, x=0 and x=Lbox are the same point!


%psi = 1 + 0.3*sin(2*pi*x/Lbox);
psi = sqrt(1 + 0.6*sin(2*pi*x/Lbox));


psiOut = zeros(length(psi),Nt+1);
VOut = zeros(length(psi),Nt+1);

psiOut(:,1) = psi;

% fourier space variables
fftw('planner','measure');
kx = (-N/2:N/2-1)' * (2*pi/Lbox);
kSq = fftshift(kx.^2);
clear kx


V = (abs(psi).^2 - 1);
V = -fftn(V);
V = V ./ ( kSq  + (kSq==0) );
V = ifftn(V);

VOut(:,1) = V;

t = 0;
dt = tFinal / Nt;
assert( dt < (1/6)*dx^2 );



%% simulation
tic;
for i=1:Nt    
    
    % kinetic - drift
    psi = fftn(psi);
    psi = exp(dt/2 * -1.i * kSq) .*psi;
    psi = ifftn(psi);
    
    % potential - kick
    V = (abs(psi).^2 - 1);
    V = -fftn(V);
    V = V ./ ( kSq  + (kSq==0) );
    V = ifftn(V);
    psi = exp(-1.i * dt * V).*psi;
    t = t + dt
    
    psiOut(:,i+1) = psi;    
    VOut(:,i+1) = V; 
    
end
toc;



end

