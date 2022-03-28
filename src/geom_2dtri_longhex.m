% function geom = geom_2dtri_longhex(K0, K1, N, alpha, bc)
%
% K0 : the defect is assumed to lie on -K0:K0
% K1 : number of atomistic layers surrounding K0
% N  : approximate radius of the domain, N - K0 - K1 is the number of
%      atomic layers surrounding the atomistic region
%      choose N = K1 for a fully atomistic simulation
% alpha : mesh growth factor
% bc : type of boundary condition, default is 'dir'
%      (if 'per' then the periodicity information is generated)
%

function geom = geom_2dtri_longhex(K0, K1, N, alpha, bc, ref_params)

if nargin == 0
  test_geom_2dtri_longhex()
  return;
end

if nargin < 5, bc = 'dir'; end

% some refinement parameters
THETA = 0.5;
KFAC = 2;
KFAC_MIN = 1;

if nargin >= 6
  THETA = ref_params(1);
  KFAC = ref_params(2);
  KFAC_MIN = ref_params(3);
end

% lattice directions, auxiliary operators
a1 = [1;0];
Q6 = [cos(pi/3), -sin(pi/3); sin(pi/3), cos(pi/3)];
a2 = Q6 * a1;
a3 = Q6 * a2;
aa = [a1, a2, a3, -a1, -a2, -a3, a1, a2];

%% atomistic core
X = a1 * (-K0:K0); nX = 2*K0+1;
c = [nX, nX, 1, 1, 1, nX];
L1 = 2*K0; L2 = 0;
for j = 1:K1
  L1 = L1 + 1; o1 = ones(1,L1);
  L2 = L2 + 1; o2 = ones(1,L2);
  
  X = [ X, ...
    (X(:,c(1)) + a1)*o2 + a3*(0:L2-1), ...
    (X(:,c(2)) + a2)*o1 - a1*(0:L1-1), ...
    (X(:,c(3)) + a3)*o2 - a2*(0:L2-1), ...
    (X(:,c(4)) - a1)*o2 - a3*(0:L2-1), ...
    (X(:,c(5)) - a2)*o1 + a1*(0:L1-1), ...
    (X(:,c(6)) - a3)*o2 + a2*(0:L2-1) ];
  
  c(1) = c(6) + L2-1; if j==1, c(1) = c(1)+1; end
  c(2) = c(1) + L2;
  c(3) = c(2) + L1;
  c(4) = c(3) + L2;
  c(5) = c(4) + L2;
  c(6) = c(5) + L1;
end

% store the indices for the core
nA = size(X, 2);
ca = c;

% number of atomic spacings left over to fill up to the boundary:
Nc = N - K0 - K1;
K = K0 + K1;

if Nc > 0

  %% create mesh spacing on one side of the hexagon
  % first create a trial spacing
  R = K;
  h_old = 1;
  while R(end) < N
    % mult = 1 - 0.5*alpha/K * (R(end)/K)^(alpha-1)
    % h = 0.5 * (R(end)/K)^alpha / mult
    h = (R(end)/K)^alpha;
    for j = 1:10
      h = THETA * (R(end)/K)^alpha + (1-THETA) * ((R(end)+h)/K)^alpha;
      % h = ((R(end)+h)/K)^alpha ;
    end
%     h = max(1, min(h, KFAC * h_old));
    h = max(KFAC_MIN * h_old, min(h, KFAC * h_old));
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
  c(7) = c(1);
  n = zeros(1,6);
  % loop through radial spacing
  for i = 2:(length(R)-1)
    % mesh size along the hexagon side
    h = 0.5 * (R(i+1)-R(i-1));
    
    % loop through the 6 sides
    for j = 1:6
      x0 = X(:,c(j)) + (R(i)-R(i-1)) * aa(:,j);
      x1 = X(:,c(j+1)) + (R(i)-R(i-1)) * aa(:,j+1);
      % estimate how many points will fit
      n(j) = ceil( norm(x1-x0) / h );
      % create scalar mesh along that edge (convex coordinates)
      % then remove last point to avoid duplication
      t = linspace(0, 1, n(j)+1);
      t = t(1:n(j));
      % create one set of mesh points ....
      X = [X, x0 * (1-t) + x1 * t];
    end
    
    c(6) = size(X,2) - n(6)+1;
    c(5) = c(6) - n(5);
    c(4) = c(5) - n(4);
    c(3) = c(4) - n(3);
    c(2) = c(3) - n(2);
    c(1) = c(2) - n(1);
    c(7) = c(1);
  end  

end


% assemble some mesh info
nX = size(X, 2);
iBdry = c(1):nX;


%% info for periodic boundary conditions!!!
%       _____nc(2)
%      /     \ side 1
%     /       \nc(1)
%     \       /
%      \_____/ side 6
%
if strcmp(bc, 'per')
  % complete corner indices
  nc = c;
  nc(7) = nX+1;
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
  geom.aper = [ N * (a2+a3), - N * (a2+a3), ...
    (N+2*K0)*a1 + N*a2, -(N+2*K0)*a1 - N*a2, ...
    (N+2*K0)*a1 - (N)*a3, -(N+2*K0)*a1 + (N)*a3 ];
  % geom.aper = N * [ a1+a2, a2+a3, a3-a1, -a1-a2, -a2-a3, -a3+a1];
end

%% create a delaunay triangulation
% first modify the domain a bit so that it becomes numerically convex
% (otherwise round-off errors cause some problems)
Y = X;
r = sqrt(sum(Y.^2, 1));
Y(:, iBdry) = Y(:, iBdry) .* ( ones(2,1) * (1 + 1e-5*(N-r(iBdry))) );
DT = DelaunayTri(Y');


onL = false * ones(1, nX);
onL(1:nA) = true;

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
geom.K0 = K0;
geom.K1 = K1;
geom.alpha = alpha;
geom.dDim = 2;
geom.id = '2dtri_longhex';
geom.plot = @geom_plot_2d;
geom.iBdry = iBdry;
geom.bc = bc;
geom.aa = aa;

% if N == K1
%   geom.fulla = true;
% else
%   geom.fulla = false;
% end
if sum(onL) == nX;
    geom.fulla = true;
else
    geom.fulla = false;
end
end



%% Test routine
function test_geom_2dtri_longhex()

% geom = geom_2dtri_longhex(3, 4, 20, 1.5, 'dir');
% geom = geom_create_vacancies(geom, 1:7);
% geom.plot(geom);

geom = geom_2dtri_longhex(3, 10, 1000, 1.5, 'dir', [0.5, 10, 1]);
% geom = geom_2dtri_longhex(3, 4, 20, 1.5, 'dir');
% geom = geom_create_vacancies(geom, 1:7);
geom = geom_analyze(geom, 2);
geom.plot(geom);
hold on;
J = find(geom.di <= 3);
plot(geom.X(1,J), geom.X(2,J), 'k*');
hold off;


% % test the periodicity information: 
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