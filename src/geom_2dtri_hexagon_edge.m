
%% function geom = geom_2dtri_hexagon_edge(N, K, alpha {, bc})
%
% generates the geometry structure for a triangular lattice, in a
% hexagonal domain, with atomistic region in the centre, EDGE DISLOCATION
%
% Input
%      N : sidelength of the hexagonal domain in atomic units
%      K : sidelength of the fully refined region (not necessarily the
%          atomistic region!)
%  alpha : mesh grading parameter: the radial mesh size function will
%          satisfy h(r) \approx (r/K)^alpha
%     bc : 'dir' or 'per'
%
% Output
%   geom : structure containing geom.N, geom.K, geom.alpha (the inputs),
%            geom.X   : vertices in A Z^2 space)
%            geom.onL : boolean array to determine whether a vertex lies
%                       on a lattice site (in the fully refined region)
%            geom.DT  : output of DelaunaTri applied to A*X
%            geom.A   : deformation matrix mapping Z^2 to the tri-lattice
%            geom.id  : '2dtri_hexagon'
%            geom.dDim  : 2  (domain dimension)
%            geom.plot  : plotting routine
%            geom.iBdry : boundary indices
%            geom.iPer  : 2 x nPer array, info for periodic b.c.s
%            geom.bc    : input (default is 'dir')
%
% WARNINIG: this routine does a lot of array-growing in loops and could
%           be slow for large problems

function geom = geom_2dtri_hexagon_edge(N, K, alpha, bc)

if nargin == 0
  test_geom_2dtri_hexagon_edge()
  return;
end

if nargin < 4
  bc = 'dir';
end

% some refinement parameters
THETA = 0.5;
KFAC = 2;

% lattice directions, auxiliary operators
a1 = [1;0];
Q6 = [cos(pi/3), -sin(pi/3); sin(pi/3), cos(pi/3)];
a2 = Q6 * a1;
a3 = Q6 * a2;
% aa = [a1, a2, a3, -a1, -a2, -a3, a1, a2];

%% atomistic core
% create a hexagon with 2*K+1 atoms across
X = zeros(2, 1 + 6*sum(1:K));
X(:,1) = [0;0]; ind = 1;
for j = 1:K
  x = Q6' * (j * a1 * ones(1, j) + a3 * (0:(j-1)));
  for k = 1:6
    x = Q6 * x;
    X(:, ind+1:ind+j) = x;
    ind = ind + j;
  end
end
nA = ind;

% in case N == K then we need to know the number of atoms on
% one sidelength of the hexagon
n = length(x);

% X = X + 5e-1 * rand(2,length(X));
X = get_2dtri_edge(X);

% % the outermost layer gets edge constraints (for constrained
% % Delaunay triangulation - this guarantees that the triangulation won't
% % destroy the lattice topology in the atomistic region.)
% iB = (nA-6*K+1):nA;
% Ca = [ [iB(1:end-1); iB(2:end)], [iB(end); iB(1)] ];


if K ~= N

  %% mesh spacing on one side of the hexagon
  % first create a trial spacing
  R = K+0.25;
  h_old = 1;
  while R(end) < N
    % mult = 1 - 0.5*alpha/K * (R(end)/K)^(alpha-1)
    % h = 0.5 * (R(end)/K)^alpha / mult
    h = (R(end)/K)^alpha;
    for j = 1:10
      h = THETA * (R(end)/K)^alpha + (1-THETA) * ((R(end)+h)/K)^alpha;
      % h = ((R(end)+h)/K)^alpha ;
    end
    h = max(1, min(h, KFAC * h_old));
    % h = max(1, h);
    R = [R, R(end) + h];
    h_old = h;
  end
  % the last point is now outside the domain
  % if R(end-1) is close to N, remove the last point and scale up
  % if R(end) is close to N, leave the last point and scale down
  h = R(end) - N;
  if h > 0.5 * (N/K)^alpha
    R = R(1:end-1);
    R = K + (R - K) * (N-K)/(R(end)-K);
  else
    R(2:end) = R(2) + (R(2:end)-R(2)) * (N-R(2))/(R(end)-R(2));
  end
  % add an extra point at the end which is needed in the next step
  R = [R, R(end) + (N/K)^alpha];
  % after rescaling some weird things can happen near K, which are
  % easily fixed as follows:
  if ( R(3)-R(2) < R(2)-R(1) )
    R(2) = 0.5 * (R(1) + R(3));
  end
  
  %% create nodes for the entire domain
  % loop through radial spacing
  for i = 2:(length(R)-1)
    % mesh size along the hexagon side
    h = 0.5 * (R(i+1)-R(i-1));
    % estimate how many points will fit
    n = ceil( R(i) / h );
    % create scalar mesh along that edge (convex coordinates)
    % (remove last point to avoid duplication)
    t = linspace(0, 1, n+1);
    t = t(1:n);
    % create one set of mesh points ....
    x = R(i) * a1 * (1-t) + R(i) * a2 * t;
    X = [X, x];
    % rotate 5 times
    for j = 1:5
      x = Q6 * x;
      X = [X, x];
    end
  end

end

% assemble some mesh info
nX = size(X, 2);
iBdry = (nX-6*n+1):nX;


%% info for periodic boundary conditions!!!
%       _____nc(2)
%      /     \ side 1
%     /       \nc(1)
%     \       /
%      \_____/ side 6
%
if strcmp(bc, 'per')
  % corner indices
  nc = (nX-6*n+1):n:(nX+1);
  % connected vertices from side 4 to side 1
  iPer = [ nc(4):nc(5); nc(2):-1:nc(1) ];
  % connected vertices from side 5 to side 2
  % we leave out nc(5)~nc(3) since nc(3)~nc(1) at the end,
  % which is already established in the first set of connections
  iPer = [iPer, ...
    [ (nc(5)+1):nc(6); (nc(3)-1):-1:nc(2) ]  ];
  % connected vertices from side 6 to side 3
  % we leave out nc(6) since nc(6)~nc(4)~nc(2), and this is already
  % established in the second set of connections
  iPer = [iPer, ...
    [ (nc(6)+1):(nc(7)-1); (nc(4)-1):-1:(nc(3)+1) ] ];
  % finally connect nc(3) with nc(1)
  iPer = [ iPer, [nc(3); nc(1)] ];
  
  % write into geom field
  geom.iPer = iPer;
  
  % periodicity vectors ("super lattice vector")
  geom.aper = N * [ a1+a2, a2+a3, a3-a1, -a1-a2, -a2-a3, -a3+a1];
end

%% create a delaunay triangulation
% first modify the domain a bit so that it becomes numerically convex
% (otherwise round-off errors cause some problems)
Y = X;
r = sqrt(sum(Y.^2, 1));
Y(:, iBdry) = Y(:, iBdry) .* ( ones(2,1) * (1 + 1e-5*(N-r(iBdry))) );
DT = DelaunayTri(Y');

% % quick test: (DEBUGGING ONLY!)
% % loop over all elements
% for k = 1:size(DT.Triangulation, 1)
%   t = DT.Triangulation(k, :);
%   J = zeros(2, 2);
%   for j = 1:2
%     J(:,j) = X(:,t(j+1)) - X(:,t(1));
%   end
%   if det(J) < 0.01
%     error('ERROR: small element in geom_2dtri_hexagon');
%   end
% end

%--------------------------------------------------------------------------
% COMMENT: should the reference configuration be Z^d or A Z^d?
%          for the moment I am taking A Z^d, since that allows us
%          to use all the auxiliary functions provided by DelaunayTri
%          but if we want Z^2, then we should uncomment the following lines
% % convert to Z^2
% X = A \ X;
% % round the lattice vertices to exact lattice sites
% X(:, onL==true) = round(X(:, onL==true));
%
% NOTE:
% for GQC it seems that we should use Z^2 as the reference configuration
% because of round-off errors in the computation of the GQC-parameters
%--------------------------------------------------------------------------


% define flags that certain nodes need to be on exact lattice sites!
nX = size(X, 2);
onL = false * ones(1, nX);
onL(1:nA-1) = true;

% write everything into the geom structure
% NOTE : X, T are copied since that allows us to make vacancies later
%        we keep DT anyways so that we can use its auxiliary functions
geom.X = X;
geom.nX = size(geom.X, 2);
geom.T = DT.Triangulation';
geom.nT = size(geom.T, 2);
geom.onL = onL;
% geom.DT = DT;
geom.A = [a1, a2];
geom.N = N;
geom.K = K;
geom.alpha = alpha;
geom.dDim = 2;
geom.id = '2dtri_hexagon_edge';
geom.plot = @geom_plot_2d;
geom.iBdry = iBdry;
geom.bc = bc;

if N == K
  geom.fulla = true;
else
  geom.fulla = false;
end

end


function test_geom_2dtri_hexagon_edge()
close all
N = 30;
K = 8;
geom = geom_2dtri_hexagon_edge(N, K, 2, 'dir');
% geom = geom_create_vacancies(geom, [1, 2]);
% geom = geom_create_vacancies(geom, [1, 2]);
% geom = geom_analyze(geom);
geom.plot(geom);
% geom.nX
% K_sq = 6*K*(K-1);
% geom.nX / K_sq


% for j = 1:length(geom.iBdry)
%   n = geom.iBdry(j);
%   text(geom.X(1,n)+0.2, geom.X(2,n), num2str(n));
% end
% for j = 1:size(geom.iPer, 2)
%   n = geom.iPer(1,j);
%   m = geom.iPer(2,j);
%   text(geom.X(1,n), geom.X(2,n)-0.5, num2str(m));
% end

end

