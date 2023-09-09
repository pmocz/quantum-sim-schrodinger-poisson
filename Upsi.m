function U = Upsi( lambda )
%UPSI Summary of this function goes here
%   Detailed explanation goes here

% n even
% depth n-1
% lambda length = 2*n*(n-1)

% https://quantaggle.com/algorithms/ansatz/
% --|-----|-----|--
%   |     |     |
% --|--|--|--|--|--
%      |     |
% --|--|--|--|--|--
%   |     |     |
% --|--|--|--|--|--
%      |     |
% --|--|--|--|--|--
%   |     |     |
% --|-----|-----|--

n = round( (1+sqrt(1+2*length(lambda)))/2 );% length(lambda);
depth = n-1;

Sx = [0 1; 1 0];
Sy = [0 -1.i; 1.i 0];
Sz = [1 0; 0 -1];
I = eye(2);
cnot = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];

% https://www.quantum-inspire.com/kbase/rotation-operators/
Rx = @(lambda) cos(lambda/2) * I + 1.i * sin(lambda/2) * Sx;
Ry = @(lambda) cos(lambda/2) * I + 1.i * sin(lambda/2) * Sy;
Rz = @(lambda) cos(lambda/2) * I + 1.i * sin(lambda/2) * Sz;


colY = cell(depth,1);
colZ = cell(depth,1);

% rotation layers
cc = 1;
for d = 1:depth
    col = Ry(lambda(cc));
    cc = cc + 1;
    colY{d} = col;
    for i = 2:n
        col = Ry(lambda(cc));
        cc = cc + 1;
        colY{d} = kron(colY{d}, col);
    end
end

for d = 1:depth
    col = Rz(lambda(cc));
    cc = cc + 1;
    colZ{d} = col;
    for i = 2:n
        col = Rz(lambda(cc));
        cc = cc + 1;
        colZ{d} = kron(colZ{d}, col);
    end
end

isOdd = (mod(n,2) == 1);

% cnot layers:
cnotA = cnot;
for i = 2:floor(n/2)
    cnotA = kron(cnotA, cnot);
end
if isOdd
    cnotA = kron(cnotA, I);
end

cnotB = I;
for i = 1:ceil(n/2)-1
    cnotB = kron(cnotB, cnot);
end
if ~isOdd
    cnotB = kron(cnotB, I);
end


% compose all the layers
U = colY{1} * colZ{1};
for d = 2:depth
    U = colY{d} * colZ{d} * cnotA * cnotB * U;
end


end

