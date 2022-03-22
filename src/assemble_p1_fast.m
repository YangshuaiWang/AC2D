
% K = assemble_p1_fast(P, T, A, c)
% Basic implementation of P1 matrix assembly, with triplet format
% to speed-up sparse matrix assembly; dimension independent.
% Input:
%   P : nP x dim array, node positions
%   T : nT x dim+1 array, mesh connectivity
%   A : kxdimx mxdim or dim x dim or scalar, diffusion tensor {default 1}
%   c : k x m or 1 x 1, coefficient for the mass matrix {default 0}
% Output:
%   K : stiffness(A) + mass(c)
%
% Note: if A and c are not symmetric, then things might get weird...

function K = assemble_p1_fast(P, T, A, c)

if nargin == 0
  test_assemble_p1_fast();
  return;
end

%% default parameters
if nargin < 3, A = 1; end

%% Some Pre-computation
% some useful constants
nP = size(P, 1); 
nT = size(T, 1); 
n = size(P, 2); 
% check dimensions
if size(T, 2) ~= n+1
  error('assemble_p1_fast : dimensions of P and T do not match.');
end
% number of components of u
m = size(c, 2);
% numebr of components in target space
k = size(c, 1);
% if the tensor A is a scalar then it is turned into the
% isotropic laplacian tensor (vectorial if necessary)
if numel(A) == 1
  if m ~= k, error('If A is scalar then c must be square.'); end
  a = A;
  A = zeros([k, n, m, n]);
  for ik = 1:k
    for in = 1:n
      A(ik,in,ik,in) = a;
    end
  end
end
% check dimension of A
if numel(A) ~= k*n*m*n, error('numel(A) must be k*n*m*n');end
A = reshape(A, [k,n,m,n]);

% some local matrices
Mloc = (eye(n+1) + ones(n+1)) / factorial(n+2);
A = A / factorial(n);

%% allocate space for matrix triplet format
rowK = zeros(k*m*(n+1)^2 * nT, 1);
colK = zeros(k*m*(n+1)^2 * nT, 1);
zK = zeros(k*m*(n+1)^2 * nT, 1);

%% assemble the matrices in triplet format
iTrip = 0;
for iT = 1:nT  
  % compute gradients of basis functions (in C) and 
  % det of the coordinate transform.
  M = [ones(n+1, 1), P(T(iT, :), :)];
  C = inv(M); C = C(2:(n+1), :);
  detJ = det(M);
  
  % assembly of the local matrix
  % loop over corners
  for r = 1:(n+1)       % col
    for s = 1:(n+1)     % row
      % loop over components
      for i = 1:m       % col
        for j = 1:k     % row
          % compute matrix coordinates of this DOF
          iTrip = iTrip + 1;

          % non-interlaced version
          % rowK(iTrip) = T(iT, s) + (j-1) * nP;
          % colK(iTrip) = T(iT, r) + (i-1) * nP;
          
          % interlaced version
          rowK(iTrip) = (T(iT,s)-1) * k + j;
          colK(iTrip) = (T(iT,r)-1) * m + i;
          
          % compute matrix entry
          Z = 0;
          for a = 1:n
            for b = 1:n
             Z = Z + A(i,a,j,b) * C(a, r) * C(b, s) + Mloc(r, s) * c(j,i);
            end
          end
          zK(iTrip) = abs(detJ) * Z;
        end
      end
    end
  end
  
end

% assemble in CCS format
K = sparse(rowK, colK, zK, k*nP, m*nP);
end


%% TEST ROUTINE
function test_assemble_p1_fast()

% try to solve laplace's equation on a uniform mesh

% mesh generation
N = 20;
x = linspace(0, 1, N+1);
[X, Y] = meshgrid(x, x);
P = [X(:), Y(:)];
iDir = find((P(:,1) == 0) | (P(:,1) == 1) | (P(:,2) == 0) | (P(:,2) == 1));
T = delaunay(P(:,1), P(:,2));

% assemble stiffness matrix
K = assemble_p1_fast(P, T, 1, 0);
K(iDir, :) = 0;
K(iDir, iDir) = speye(length(iDir));

% assemble mass matrix
M = assemble_p1_fast(P, T, 0, 1);

% assemble rhs
F = M * ones((N+1)^2, 1);
F(iDir) = 0;

% solve linear system
u = K \ F;
trisurf(T, P(:,1), P(:,2), u);

end



