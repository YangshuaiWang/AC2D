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

function geom = geom_2dtri_frenkel(K1, N, rCutH, alpha, ref_params)

if nargin == 0
  test_geom_2dtri_frenkel()
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
geom = geom_analyze(geom, rCutH);
di = geom.di;

% move one atom to another cell
iV = 14;
xI = [1.7; 0.3];

% create a vacancy
geom = geom_create_vacancies(geom, iV);
di(iV) = [];
geom.di = di;

% create an interstitial
geom = geom_create_interstitial(geom, xI);
di = geom.di;

% some more data
geom = geom_analyze(geom, rCutH);
geom.K = K1;
geom.N = N;
geom.di = di;

end


%% Test routine
function test_geom_2dtri_frenkel()

geom = geom_2dtri_frenkel(5, 15, 1);
geom.plot(geom);

end