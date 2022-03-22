% geom = geom_2dtri_frenkel(K1, N, rCutH, alpha)
%
% K0 : number of vacancy sites
%      if nV = 2*K0+1, then we remove {-L, ..., L}
%      if nV = 2*K0, then we remove {-K0+1,...,K0}
% K1 : number of atomistic layers surrounding the vacancies
% N  : approximate radius of the domain, N - K0 - K1 is the number of
%      atomic layers surrounding the atomistic region
%      choose N = K1 for a fully atomistic simulation
% rCutH : 
%

function geom = geom_2dtri_interstitial(K1, N, rCutH, alpha, ref_params)

if nargin == 0
  test_geom_2dtri_interstitial()
  return;
end

if nargin < 4
  alpha = 1.5;
end

if nargin < 5
  ref_params = [0.5, 3, 1.0];
end

% generate geometry
geom = geom_2dtri_longhex(0, K1, N, alpha, 'dir', ref_params);

% redo the triangulation so we have access to the nice methods
Y = geom.X;
r = sqrt(sum(Y.^2, 1));
Y(:, geom.iBdry) = Y(:, geom.iBdry) .* ( ones(2,1) * (1 + 1e-5*(N-r(geom.iBdry))) );
DT = delaunayTriangulation(Y');

% position of interstitial
xI = [0.5; 0.0];

% find the two neighbouring elements and their corners
iT_I = pointLocation(DT, [xI' + [0, 0.1]; xI' - [0, 0.1]]);

% remove the two elements
T = DT.ConnectivityList';
T(:, iT_I) = [];

% add the interstitial to the lattice
geom.X = [geom.X, xI];
geom.nX = size(geom.X, 2);
geom.onL(geom.nX) = true;
nX = geom.nX;

% construct the new element, connecting the interstitial to the rest
T = [T, [1;nX;3], [nX;2;3], [1;7;nX], [nX;7;2] ];
geom.T = T;
geom.nT = size(geom.T, 2);

% create neighbourhood information
geom = geom_analyze(geom, rCutH);
geom.K = K1;
geom.nV = 0;

end


%% Test routine
function test_geom_2dtri_interstitial()

geom = geom_2dtri_interstitial(6, 20, 1);
geom.plot(geom);

end