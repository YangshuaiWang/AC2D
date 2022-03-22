
%% function r = get_nhd_hop(A, rCutH)
%
% computes an interaction neighbourhood in the lattice A Z^d that can be 
% used, e.g., to compute the Cauchy--Born potential using Wcb
%
% Input:
%       A : dxd matrix, A * Z^d is a homogeneous deformed state
%           e.g. the ground state of a crystal, or something nearby
%   rCutH : cut-off radius of the interaction in hopping distance
%
% Output:
%       r : d x M, set of all A Z^d lattice vectors r s.t. with
%           hopping distance rCutH or less from the origin
%
% WARNING: no physical definition of nearest neighbours is used here
%          and the nn-relation is simply given by the delaunay 
%          triangulation of the lattice A Z^d                (TODO)
%          (since this is only a preliminary code, maybe we don't
%           want to worry about this)
%
% WARNING: this routine is slow and should not be used in loops
% NOTE   : in later versions, where hopping neighbourhoods are 
%          replaced with physical nhds, look at get_nhd in the
%          backup directory
%

function r = get_nhd_hop(A, rCutH)

if nargin == 0
  test_get_nhd_hop();
  return;
end

%% read input
% dimension of problem
d = size(A, 2);
% get smalles singular value of A
s0 = min(svd(A));

% maximum inf-norm of Z^d vectors that we want to get
rmax = ceil(rCutH / s0 * 1.5);

%% compute preliminary set of vectors
% compute larger Z^d region
t = (-rmax):rmax;
switch d
  case 1
    r = t;
  case 2
    [rx, ry] = meshgrid(t,t);
    r = [rx(:)'; ry(:)'];
  case 3
    [rx, ry, rz] = meshgrid(t,t,t);
    r = [rx(:)'; ry(:)'; rz(:)'];
end
% transform
r = A * r;
% find origin
I0 = find( (~r(1,:)) & (~r(2,:)) );

%% create neighbour matrix
% get delaunay triangularion
DT = DelaunayTri(r(1,:)', r(2,:)');
T = DT.Triangulation;
% nearest-neighbour matrix
if d ~= 2, error('ERROR: get_nhd_hop only implemented for dDim = 2'); end
NN = sparse( [T(:,1); T(:,1); T(:,2); T(:,2); T(:,3); T(:,3)], ...
             [T(:,2); T(:,3); T(:,1); T(:,3); T(:,1); T(:,2)], ...
             0.1*ones(6*size(T,1), 1) );
NN = ceil(NN);
% rCutH neighbour matrix
Neigs = NN;
for j = 2:rCutH
  Neigs = Neigs * NN;
end
% remove self-interaction
for j = 1:size(NN, 1), Neigs(j,j) = 0; end

%% get neighbourhood of origin
iN = find(Neigs(:, I0));  %#ok
r = r(:, iN);             %#ok


end


%% TEST CODE
function test_get_nhd_hop()
% triangular lattice
A = [ [1;0] [cos(pi/3); sin(pi/3)] ];
% two hops
rCutH = 3;
% compute nhd
tic, r = get_nhd_hop(A, rCutH); toc
% plot
close all;
plot(r(1,:), r(2,:), 'bo');
end
