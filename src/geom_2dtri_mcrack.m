% function geom = geom_2dtri_longhex(K0, K1, N, rCutH)
%
% K0 : the defect is assumed to lie on -K0:K0
% K1 : number of atomistic layers surrounding K0
% N  : approximate radius of the domain, N - K0 - K1 is the number of
%      atomic layers surrounding the atomistic region
%      choose N = K1 for a fully atomistic simulation
% rCutH : 
%

function geom = geom_2dtri_mcrack(K0, K1, N, rCutH, alpha)

if nargin == 0
  test_geom_2dtri_mcrack()
  return;
end

if nargin < 5
  alpha = 1.5;
end

% generate geometry
geom = geom_2dtri_longhex(K0, K1, N, alpha, 'dir');

% remove the crack
geom = geom_create_vacancies(geom, 1:(2*K0+1));

% create neighbourhood information
geom = geom_analyze(geom, rCutH);


end


%% Test routine
function test_geom_2dtri_mcrack()

geom = geom_2dtri_mcrack(3, 4, 20, 2);
geom.plot(geom);

end