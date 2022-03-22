% P = prec_laplace(geom, model{, aa})
%
% INPUT
%    geom, model : MAC geometry and model structures; uses 
%                  geom.X, P, nX, iBdry, iPer
%                  model.bc, rDim
%             aa : 2x1 array, P = aa(1)*diam * K + aa(2)/diam * M
%                  default: aa = [0.3, 1e-3];
% OUPUT
%   P : a1 * diam * K + a2/diam * M where K is the stiffness and M the 
%       mass matrix
%

function P = prec_laplace(geom, model, aa)

if nargin == 0
  test_prec_laplace();
  return;
end

% default parameters
if nargin < 3
  aa = [0.3, 1e-3];
end

% domain diameter for rescaling
d = max( max(geom.X(1,:)) - min(geom.X(1,:)), ...
         max(geom.X(2,:)) - min(geom.X(2,:)) );
% assemble stiffness matrix
P1 = assemble_p1_fast(geom.X', geom.T', aa(1) * d, aa(2)/d);

% trick apply_bc() into believing that rDim = 1 and apply the bc
switch model.bc
  case ''
    % neumann >> do nothing!
    iFree = 1:geom.nX;
  case 'dir'
    if isfield(geom, 'iDir')
      iFree = setdiff(1:geom.nX, geom.iDir);
    else
      iFree = setdiff(1:geom.nX, geom.iBdry);
    end
    P1 = P1(iFree, iFree);
    
  case 'per'
    for j = 1:size(geom.iPer, 2)
      P1(:,geom.iPer(2,j)) = P1(:,geom.iPer(2,j)) + P1(:,geom.iPer(1,j));
    end
    P1(:, geom.iPer(1,:)) = [];
    for j = 1:size(geom.iPer, 2)
      P1(geom.iPer(2,j),:) = P1(geom.iPer(2,j),:) + P1(geom.iPer(1,j),:);
    end
    P1(geom.iPer(1,:), :) = [];
    iFree = setdiff(1:geom.nX, geom.iPer(1,:));
  otherwise
    error('ERROR: prec_laplace encountered an unknown bc');
end

% combine copies of P1 for vectorial problems
nFree = length(iFree);
P = spalloc(model.rDim*nFree, model.rDim*nFree, model.rDim*nnz(P1));
for j = 1:model.rDim
  inds = j:model.rDim:(model.rDim*nFree);
  P(inds, inds) = P1;  %#ok
end

end






%% TEST FUNCTION
function test_prec_laplace()
% generate model and geometry
model = model_morse(4, 1, 2);
geom = geom_2dtri_hexagon(20, 3, 1.1, 'per');
geom = geom_analyze(geom);
model.rDim = 1;

% mass matrix
model.bc = '';
M = prec_laplace(geom, model, [0, 1]);

% Problem 1: Dirichlet
model.bc = 'dir';
A = prec_laplace(geom, model, [1, 0]);
F = M * ones(geom.nX, 1);
F(geom.iBdry) = [];
U = zeros(geom.nX, 1);
iFree = setdiff(1:geom.nX, geom.iBdry);
U(iFree) = A \ F;
figure(1);
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', U);
title('Dirichlet');

% Problem 2: Periodic
model.bc = 'per';
B = prec_laplace(geom, model, [1, 0.1]);
r2 = (geom.X(1,:).^2 + geom.X(2,:).^2)';
G = M * exp(- 5 * r2 / geom.N);
G(1) = 5; G(10) = -5;
for j = 1:size(geom.iPer, 2)
  G(geom.iPer(2,j)) = G(geom.iPer(2,j)) + G(geom.iPer(1,j));
end
G(geom.iPer(1,:)) = [];

V = zeros(geom.nX, 1);
iFree = setdiff(1:geom.nX, geom.iPer(1,:));
V(iFree) = B \ G;
V(geom.iPer(1,:)) = V(geom.iPer(2,:));
figure(2);
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', V);
end
