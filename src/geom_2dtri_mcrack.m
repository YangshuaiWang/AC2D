% function geom = geom_2dtri_longhex(K0, K1, N, rCutH)
%
% K0 : the defect is assumed to lie on -K0:K0
% K1 : number of atomistic layers surrounding K0
% N  : approximate radius of the domain, N - K0 - K1 is the number of
%      atomic layers surrounding the atomistic region
%      choose N = K1 for a fully atomistic simulation
% rCutH : 
%

function geom = geom_2dtri_mcrack(nV, K1, N, rCutH, alpha, ref_params)

if nargin == 0
  test_geom_2dtri_mcrack()
  return;
end

if nargin < 5
  alpha = 1.5;
end

if nargin < 6
  ref_params = [0.5, 3, 1.0];
end

% CASE 1: odd number of vacancy sites
if mod(nV, 2) == 1
  K0 = (nV-1)/2;
else
  % CASE 2: even number of vacancy sites
  K0 = nV/2;
end

% % generate geometry
% geom = geom_2dtri_longhex(K0, K1, N, alpha, 'dir');
% generate geometry
geom = geom_2dtri_longhex(K0, K1, N, alpha, 'dir', ref_params);
geom = geom_analyze(geom, rCutH);
di = geom.di;

% % remove the crack
% geom = geom_create_vacancies(geom, 1:(2*K0+1));
% remove the crack
if nV > 0
  di(1:nV) = [];
  geom = geom_create_vacancies(geom, 1:nV);
end

% create neighbourhood information
geom = geom_analyze(geom, rCutH);
geom.di = di;
geom.K = K1;
geom.nV = nV;

end


%% Test routine
function test_geom_2dtri_mcrack()

geom = geom_2dtri_mcrack(3, 4, 20, 2);
geom.plot(geom);

end