
% function geom = bqc_prep_geom(geom, max_dh, blendw, blendfcn)
%
% prepares geom so that it can be used to run bqc_energy
%
% Input
%     geom : valid geometry structure
%   max_dh : maximal hopping distance for atomistic interaction
%       bw : blending width in site hops
%  blendfc : string, type of blending function. Current choices:
%            'linearH', 'cubicH'
%
% Output
%   geom : * Added fields max_dh, blendw, which are simply the values given
%          * Effective volumes of elements are stored in geom.volT
%          * Effective volumes of vertices: 1 in true atomistic region,
%            zero in true continuum region, blended inbetween;
%            stored in geom.volX
%


function geom = bqc_prep_geom(geom, max_dh, blendw, blendfcn)

if nargin == 0
  test_bqc_prep_geom();
  return;
end

% check whether geom is valid > we need that geom_analyze has already
% been run on this
if (~isfield(geom, 'geom_analyze') || ~geom.geom_analyze )
  error('ERROR: must call geom_analyze before bqc_prep_geom!');
end 

%% 1. easy part
if (max_dh < 1) 
  error('ERROR: bqc_prep_geom requires max_dh >= 1');
end
if ~strcmp(blendfcn, 'qce') && blendw < 1
  error('ERROR: bqc_prep_geom requires blendw >= 1, unless qce is used');
end
geom.max_dh = max_dh;
geom.blendw = blendw;

%% 2. Compute Effective volumes of vertices
% set all volX = 0, then set them to the right value in the core atomistic 
% region
geom.volX = zeros(1, geom.nX);
% call the blending functions
Spl1 = @(s)(s);
Spl3 = @(s)(3 * s.^2 - 2 * s.^3);
switch blendfcn
  case 'cb'
    % do nothing
  case 'qce'
    geom = blend_qce(geom);
  case 'linearH'
    geom = blend_spline_hop(geom, Spl1);
  case 'cubicH'
    geom = blend_spline_hop(geom, Spl3);
  case 'laplace'
    geom = blend_laplace(geom, Spl1);
  case 'laplace+cubic'
    geom = blend_laplace(geom, Spl3);
  case 'phys-cubic'
    geom = blend_phys_spline(geom, Spl3);
  case 'phys-linear'
    geom = blend_phys_spline(geom, Spl1);
  case 'min-hessian'
    geom = blend_min_hessian(geom, false);
  case 'min-hessian-bench'
    % if we want the benchmark version, so the linear blending first
    % and use it to compute the blending region
    geom = blend_phys_spline(geom, Spl1);
    % true means, that it is a benchmark run
    geom = blend_min_hessian(geom, true);
  case 'min-hessian-bench-new'
    % if we want the benchmark version, so the linear blending first
    % and use it to compute the blending region
    geom = blend_phys_spline(geom, Spl1);
    % true means, that it is a benchmark run
    geom = blend_min_hessian_new(geom, true);
%   case 'bilaplace'
%     geom = blend_bilaplace(geom, false);
  case 'gf-min-L2'
    error('ERROR : bqc_prep_geom, gf_min_L2 blending not implemented');
  case 'gf-min-H-1'
    error('ERROR : bqc_prep_geom, gf_min_H-1 blending not implemented');
  case 'hessian'
    error('ERROR: bqc_prep_geom - hessian blending is not implemented');
  otherwise
    error('ERROR: unkown blending fcn in bqc_prep_geom');
end


%% 3. compute volT
% the second piece of information we need are the effective volumes
% of the elements T. We define
%     v_T = \sum_{a \in L} volX(a) * vol( V_a \cap T ),
% where V_a is the voronoi cell associated with the lattice site a

% write this into geom.volT(k), k = 1...nT
geom.volT = zeros(1, geom.nT);
facd = factorial(geom.dDim);
detA = det(geom.A);
for k = 1:geom.nT
  
  % TODO: calculate physical volume of T
  p = geom.X(:, geom.T(:, k));
  J = zeros(geom.dDim, geom.dDim);
  for j = 1:geom.dDim
    J(:, j) = p(:, j+1) - p(:, 1);
  end
  geom.volT(k) = abs(det(J)) / facd;
  
  % in Z^d each Voronoi cell has volume = 1, it has volume = det(A)
  % in A Z^d. But in the QCE energy, we want to keep the effective
  % volume associated with each lattice site at 1, hence we need to
  % rescale the element volumes:
  geom.volT(k) = geom.volT(k) / detA;
  
  % add up the effective volumes of the vertices of T
  % >> get nA, nAext
  vA = 0; nAext = 0;
  for j = 1:(geom.dDim+1)
    n = geom.T(j, k);
    vA = vA + geom.volX(n) / 6;
    if geom.di(n) >= 0
      nAext = nAext + 1;
    end
  end
  
  % TODO: the above loop makes a huge assumption!!! namely, that the
  %       volume of each Voronoi cell within T is the same!!!!
  %       we need to fix this part and compute the effective volumes 
  %       properly. The current version probably works only in the 2D
  %       triangular lattice.

  
  % if there is at least one such vertex, then we require in the current
  % version of the code that the element T must be a nearest-neighbour
  % element, i.e. all corners should belong to the extended atomistic
  % region.
  if ( (vA > 0) && (nAext ~= geom.dDim+1) )
    error(['ERROR: bqc_prep_geom has encountered a situation ' ...
           'it cannot cope with yet']);
  end
  
  % adjust volume of the element, and correct for round-off errors
  geom.volT(k) = max(0, geom.volT(k) - vA);
  if geom.volT(k) < 1e-12
    geom.volT(k) = 0;
  end
 
end

end


%% function geom = blend_qce(geom)
function geom = blend_qce(geom)
geom.volX(geom.di >= geom.max_dh) = 1;
end

%% function geom = blend_spline_hop(geom, Spl)
% generates a blending function that is a polynomial function of the
% hopping distance.
function geom = blend_spline_hop(geom, Spl)
% loop through vertices
for n = 1:geom.nX
  if geom.di(n) >= (geom.max_dh + geom.blendw)
    geom.volX(n) = 1;
  elseif geom.di(n) >= geom.max_dh
    % evaluate cubic spline
    geom.volX(n) = Spl( (geom.di(n) - geom.max_dh + 1) / geom.blendw );
  end
end
end


%% function geom = blend_laplace(geom)
% generates a blending function by minimizing the laplacian functional
% subject to boundary condition b = 1 in the atomistic region and
% b = 0 in the continuum region
function geom = blend_laplace(geom, Spl)
K = assemble_p1_fast(geom.X', geom.T', 1, 0);
F = zeros(geom.nX, 1);
% Dirichlet nodes, where we set b = 1
iOne = find(geom.di >= geom.max_dh + geom.blendw - 1);
K(iOne, :) = 0; 
K(iOne, iOne) = speye(length(iOne));
F(iOne) = 1;
% Dirichlet nodes where we set b = 0
iZero = find(geom.di <= geom.max_dh-1);
K(iZero, :) = 0;
K(iZero, iZero) = speye(length(iZero));

% solve for blending function
geom.volX = (K \ F)';

% now smooth them using the spline Spl
for n = 1:geom.nX
  geom.volX(n) = Spl(geom.volX(n));
end

end


%% function geom = blend_phys_spline(geom)
% generates a blending function that is a polynomial function of the
% of a physical distance from the origin
% WARNING : this makes only sense for atomistic regions that
%           are radial about the origin
function geom = blend_phys_spline(geom, Spl)
if ~strcmp(geom.id, '2dtri_hexagon')
  error(['ERROR: blend_phys_spline should only be used for regions' ...
    ' created using geom_2dtri_hexagon.m in the gf-benchmark']);
end
s_inn = (geom.K - 2 - geom.blendw) * sqrt(3)/2;
s_out = (geom.K - 2) * sqrt(3)/2;
s = sqrt(sum(geom.X.^2, 1));
s = 1 - (s - s_inn) / (s_out - s_inn);
geom.volX = Spl(s); % 3 * s.^2 - 2 * s.^3;
geom.volX(s > 1) = 1;
geom.volX(s < 0) = 0;
end


%% function geom = blend_min_hessian(geom)
% computes a blending function by minimizing the L2-norm of a 
% discrete hessian
% WARNING: this implementation assumes that the mesh is fully
%          refined in the blending region
%
function geom = blend_min_hessian(geom, gfbench)
% atomistic, blending nodes, and continuum
if gfbench
  % GF-Benchmark Case (this is only used to give a fair comparison
  % against the radial cubic and linear splines
  I1 = find(geom.volX == 1);
  I0 = find(geom.volX == 0);
  iB = find(geom.di >= 1);    % this is inefficient, but only a benchmark
else
  % General case (this will normally be used)
  mx_dh = max(2, geom.max_dh);
  I0 = find(geom.di < mx_dh);
  I1 = find(geom.di >= mx_dh + geom.blendw-1);
  iB = find( (geom.di < mx_dh + geom.blendw) ...
           & (geom.di >= mx_dh - 1) );
end

% compute matrix to be minimized
H1 = pseudo_laplace(geom, iB);


% We want to minimize || H b ||^2, so we actually need H' * H
H = H1' * H1;
% finally, we want to apply boundary conditions
Ibc = [I1, I0];
H(Ibc, :) = 0; H(Ibc, Ibc) = speye(length(Ibc));
% right-hand side for boundary conditions
rhs = zeros(geom.nX, 1);
rhs(I1) = 1;
% solve linear system
geom.volX = ( H \ rhs )';

if ( (max(geom.volX) > 1) || (min(geom.volX) < 0) )
  error('bqc_prep_geom : blend_min_hessian : volX \notin [0, 1] !!');
end

% %% debugging
% b1 = geom.volX;
% Spl = @(s)(3 * s.^2 - 2 * s.^3);
% b2 = Spl(b1);
% b2(b1 == 1) = 1; b2(b1 == 0) = 0;
% subplot(2,2,1);
% trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', b1');
% title(['||Hb1|| ~ ', num2str(norm(H1*b1', 2))]);
% N = geom.N;
% axis([-N, N, -N, N, -1, 2]);
% subplot(2,2,2);
% trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', b2');
% title(['||Hb2|| ~ ', num2str(norm(H1*b2', 2))]);
% axis([-N, N, -N, N, -1, 2]);
% pause;
end


%% function geom = blend_min_hessian_new(geom)
% computes a blending function by minimizing the L2-norm of a 
% discrete hessian --- Brian's version
% WARNING: this implementation assumes that the mesh is fully
%          refined in the blending region
%
function geom = blend_min_hessian_new(geom, gfbench)
% atomistic, blending nodes, and continuum
if gfbench
  % GF-Benchmark Case (this is only used to give a fair comparison
  % against the radial cubic and linear splines
  I1 = find(geom.volX == 1);
  I0 = find(geom.volX == 0);
  % all edges inside the blending region
  iE = find( (mod(geom.volX(geom.E(1,:)), 1) > 0) | ...
             (mod(geom.volX(geom.E(2,:)), 1) > 0) );
  % the where the pseudo-laplacian (for stabilisation)
  % is assembled at
  iB = find(geom.di >= 1);
else
  % General case (this will normally be used)
  error('ERROR: the general case for blend_min_hessian_new is not implemented yet.');
%   I0 = find(geom.di < geom.max_dh);
%   I1 = find(geom.di >= geom.max_dh + geom.blendw-1);
%   iB = find( (geom.di < geom.max_dh + geom.blendw) ...
%            & (geom.di >= geom.max_dh-1) );
end

% assemble matrix to be minimized.
% We want to minimize || H1 b ||^2, but H1'*H1 is singular, so we
% minimize || H1 b ||^2 + epsilon || H2 b ||^2 instead.
H1 = bvk_matrix(geom, iE);
H2 = pseudo_laplace(geom, iB);
H = H1' * H1 + 1e-2 * H2'*H2;

% apply boundary conditions
Ibc = [I1, I0];
H(Ibc, :) = 0; H(Ibc, Ibc) = speye(length(Ibc));
% right-hand side for boundary conditions
rhs = zeros(geom.nX, 1);
rhs(I1) = 1;
% solve linear system
geom.volX = ( H \ rhs )';

end


%% sort_neigs
% returns nearest-neighbours, ordered in anti-clockwise order.
% gives an error is there are not 6 neighbours. Probably also 
% if the six neighbours are not precisely on a triangular lattice, then
% the code will crash.
%
function iN = sorted_neigs(i0, geom)
iN = find(geom.NN(:, i0));
if length(iN) ~= 6
  iN = [];
  return;
  % error('ERROR: sorted_neigs assumes there are 6 neighbours!');
end

x = geom.X(:, iN) - geom.X(:,i0) * ones(1,6);
x1 = x(:,1); x1p = [-x1(2); x1(1)];

ii = zeros(1,6); ii(1) = 1;
for i = 2:6
  if x1'*x(:,i) > 0
    if x1p'*x(:,i) > 0
      ii(2) = i;
    else
      ii(6) = i;
    end
  else
    b = x1p'*x(:,i);
    if abs(b) < 1e-5
      ii(4) = i;
    elseif b > 0
      ii(3) = i;
    else 
      ii(5) = i;
    end
  end
end
iN = iN(ii);

end



%% H1 = pseudo_laplace(geom, iB)
% compute a matrix H1 such that H1 * B is a sort of discrete
% second derivative in each component. This is used to generate
% minimization problems of the norm    min |H1*B|^2
% subject to some constraints
function H1 = pseudo_laplace(geom, iB)
% number of triplets: 9 per vertex where a second derivative is computed
ntrip = 9 * length(iB);
% allocate space
iRow = zeros(ntrip, 1);
iCol = zeros(ntrip, 1);
Val = zeros(ntrip, 1);
% assemble second-derivative operator in triplet format
iT = 0;     % (iT is the current triplet index)
iR = 0;     % (iR is the current row-index)
mod6 = @(i_)(mod(i_-1, 6)+1);  % modified modulus function
for j = 1:length(iB)
  % get the index of the current node in geom.X
  i = iB(j);
  % get sorted neighbour list
  iN = sorted_neigs(i, geom);
  % compute three second derivative functionals and add them to H1
  if ~isempty(iN)
    for k = 1:3
      iT = iT + 1; iR = iR + 1;
      iRow(iT) = iR; iCol(iT) = i; Val(iT) = -2;
      iT = iT + 1;
      iRow(iT) = iR; iCol(iT) = iN(k); Val(iT) = 1;
      iT = iT + 1;
      iRow(iT) = iR; iCol(iT) = iN(mod6(k+3)); Val(iT) = 1;
    end
  end
end
% assemble the operator into CCS sparse matrix
H1 = sparse(iRow(1:iT), iCol(1:iT), Val(1:iT), iR, geom.nX);
end


%% H1 = bvk_matrix(geom, iE)
% each row of H1 represents a linear functional of the form
%    (H1 * B)_i = D_j D_{j+1} B_i, which is more or less a mixed
% second derivative, and comes naturally out of minimizing the ghost
% force. See note
%      2011-08-11-bvk-BQC_Ghost_F.pdf
%
function H1 = bvk_matrix(geom, iE)
% number of triplets: 9 per vertex where a second derivative is computed
ntrip = 4 * length(iE);
% allocate space
iRow = zeros(ntrip, 1);
iCol = zeros(ntrip, 1);
Val = zeros(ntrip, 1);
% assemble second-derivative operator in triplet format
iT = 0;  % (iT is the triplet index of the matrix)
iR = 0;  % (row-index of matrix)
for j = 1:length(iE)
  % nodes that get a -1 coefficient
  im = geom.E(:, iE(j));
  % nodes that get a +1 coefficient: use edge-element adjacancy relation
  % to get the vertex opposite the edges in the two neighbouring elements
  ip = [0,0];
  for r = 1:2
    t = geom.aET(r, iE(j));
    ip(r) = setdiff(geom.T(:, t), im);
  end
  % create four triplets
  iR = iR + 1;
  iT = iT+1; iRow(iT) = iR; iCol(iT) = im(1); iVal(iT) = -1;
  iT = iT+1; iRow(iT) = iR; iCol(iT) = im(2); iVal(iT) = -1;
  iT = iT+1; iRow(iT) = iR; iCol(iT) = ip(1); iVal(iT) = 1;
  iT = iT+1; iRow(iT) = iR; iCol(iT) = ip(2); iVal(iT) = 1;
end
% assemble the operator into CCS sparse matrix
H1 = sparse(iRow, iCol, Val, iR, geom.nX);
end


%% test code
function test_bqc_prep_geom()
figure(1);
N = 50; dH = 3; bw = 5; K = dH + 2*bw;
% profile on;
geom = geom_2dtri_hexagon(N, K, 3/2);
geom = geom_analyze(geom);
geom = bqc_prep_geom(geom, dH, bw, 'min-hessian');
% profile off;
% profile report;
geom.plot(geom);
% for n = 1:geom.nX
%   text(geom.X(1,n)+0.1, geom.X(2,n), num2str(geom.volX(n)));
% end
% for k = 1:geom.nT
%   % get barycentre of element k
%   x = 0.333 * sum(geom.X(:, geom.T(:, k)), 2);
%   text(x(1), x(2), num2str(geom.volT(k), '%4.2f'));
% end

figure(2);
% plot blending function
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', geom.volX);
axis([-N, N, -N, N, -2, 3]);
end











%% BACKUP CODES


% %% OLD CODE
% %% function geom = blend_min_hessian(geom)
% % computes a blending function by minimizing the L2-norm of a 
% % discrete hessian
% % WARNING: this implementation assumes that the mesh is fully
% %          refined in the blending region
% %
% function geom = blend_min_hessian(geom, gfbench)
% % atomistic, blending nodes, and continuum
% if gfbench
%   % GF-Benchmark Case (this is only used to give a fair comparison
%   % against the radial cubic and linear splines
%   I1 = find(geom.volX == 1);
%   I0 = find(geom.volX == 0);
%   iB = find(geom.di >= 1);    % this is inefficient, but only a benchmark
% else
%   % General case (this will normally be used)
%   I0 = find(geom.di < geom.max_dh);
%   I1 = find(geom.di >= geom.max_dh + geom.blendw-1);
%   iB = find( (geom.di < geom.max_dh + geom.blendw) ...
%            & (geom.di >= geom.max_dh-1) );
% end
% % number of triplets: 7 per vertex where a second derivative is computed
% ntrip = 7 * length(iB);
% % allocate space
% iRow = zeros(ntrip, 1);
% iCol = zeros(ntrip, 1);
% Val = zeros(ntrip, 1);
% % assemble second-derivative operator in triplet format
% iT = 0;  % (iT is the triplet index)
% for j = 1:length(iB)
%   i = iB(j);
%   % get neighbour list
%   iN = find(geom.NN(:, i));
%   if length(iN) ~= 6
%     error('ERROR: blend_min_hessian assumes there are 6 neighbours!');
%   end
%   % 
%   iT = iT + 1;
%   iRow(iT) = i; iCol(iT) = i; Val(iT) = -6;
%   for k = 1:length(iN)
%     iRow(iT+k) = i; iCol(iT+k) = iN(k); Val(iT+k) = 1;
%   end
%   iT = iT + 6;
% end
% % assemble the operator into CCS sparse matrix
% H1 = sparse(iRow, iCol, Val, geom.nX, geom.nX);
% % We want to minimize || H b ||^2, so we actually need H' * H
% H = H1' * H1;
% % finally, we want to apply boundary conditions
% Ibc = [I1, I0];
% H(Ibc, :) = 0; H(Ibc, Ibc) = speye(length(Ibc));
% % right-hand side for boundary conditions
% rhs = zeros(geom.nX, 1);
% rhs(I1) = 1;
% % solve linear system
% geom.volX = ( H \ rhs )';
% 
% % %% debugging
% % b1 = geom.volX;
% % Spl = @(s)(3 * s.^2 - 2 * s.^3);
% % b2 = Spl(b1);
% % b2(b1 == 1) = 1; b2(b1 == 0) = 0;
% % subplot(2,2,1);
% % trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', b1');
% % title(['||Hb1|| ~ ', num2str(norm(H1*b1', 2))]);
% % N = geom.N;
% % axis([-N, N, -N, N, -1, 2]);
% % subplot(2,2,2);
% % trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', b2');
% % title(['||Hb2|| ~ ', num2str(norm(H1*b2', 2))]);
% % axis([-N, N, -N, N, -1, 2]);
% % pause;
% end






% % PDCO SOLUTION
% after this preparation, we want to minimize
%    | H1 * b |^2 subject to the constraint that b = 0,1 in A,C resp.
% but this seems to be ill-posed (H1'*H1 is not invertible)
% so instead, we reformulate the problem:
%      Minimize   - sum_j b_j      (to make b_j as close to 1 as possible)
%    subject to      * H1'*H1 b = 0
%                    * 0 <= b <= 1
%                    * b = 0 in continuum region
%                    * b = 1 in atomistic region
% % PREPARE FOR PDCO
% addpath ./packages/pdco
% % linear constraint
% pdMat = H1' * H1;
% b = zeros(geom.nX, 1);
% % upper and lower bounds
% bl = zeros(geom.nX, 1); bl(I1) = 1;
% bu = ones(geom.nX, 1); bu(I0) = 0;
% % objective
% pdObj = -ones(geom.nX, 1);
% % regularizations
% d1 = 1e-4; d2 = 1e-4;
% % options
% pdco_opts = pdcoSet();
% % call pdco
% [x, y, z, inform] = pdco(pdObj, pdMat, b, bl, bu, d1, d2, pdco_opts, ...
%   geom.volX', zeros(geom.nX, 1), zeros(geom.nX, 1), 1, 1);




% % OLD VERSION: regularized problem - gives no useful solution
% % We want to minimize || H b ||^2, so we actually need H' * H
% H = H1' * H1 + 1e-6 * speye(geom.nX);
% % finally, we want to apply boundary conditions
% Ibc = [I1, I0];
% H(Ibc, :) = 0; H(Ibc, Ibc) = speye(length(Ibc));
% % right-hand side for boundary conditions
% rhs = ones(geom.nX, 1);
% rhs(I1) = 0;
% % solve linear system
% geom.volX = 1 - ( H \ rhs )';



% %% function geom = blend_bilaplace(geom)
% % generates a blending function by minimizing the laplacian functional
% % subject to boundary condition b = 1 in the atomistic region and
% % b = 0 in the continuum region
% function geom = blend_bilaplace(geom)
% K = assemble_p1_fast(geom.X', geom.T', 1, 0);
% F = zeros(geom.nX, 1);
% % Dirichlet nodes, where we set b = 1
% % OLD: iOne = find(geom.di >= geom.max_dh + geom.blendw - 1);
% iOne = find(geom.volX == 1);
% K(iOne, :) = 0;
% % Dirichlet nodes where we set b = 0
% % OLD: iZero = find(geom.di <= geom.max_dh-1);
% iZero = find(geom.volX == 0);
% K(iZero, :) = 0;
% % with K as above we now have that 
% %    K * b ~ - \Delta b  at the blending nodes, and zero otherwise
% % to minimize || K * b || we solve
% %    K' * K * b = 0
% % subject to boundary conditions
% K = K' * K;
% K(iOne, :) = 0; K(iOne, iOne) = speye(length(iOne));
% K(iZero, :) = 0; K(iZero, iZero) = speye(length(iZero));
% F(iOne) = 1;
% 
% % solve for blending function
% geom.volX = (K \ F)';
% 
% end


