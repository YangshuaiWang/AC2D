% function geom = geom_refneigs(geom, rCutH)
%
% generates an interaction neighbourhood in the reference configutation
% useful for static neighbourhood computations and for intial Verlet list
% setup
%

function geom = geom_refneigs(geom, rCutH)

if nargin == 0
  test_geom_refneigs();
  return;
end

% generate the complete interaction neighbourhood
% WARNING: in the current version, this is static, but in 
% future versions this must be dynamic
Ng = 0.1 * geom.NN;
for j = 1:(rCutH-1)
  Ng = Ng + geom.NN*Ng;
end
Ng = ceil(Ng);

for j = 1:size(Ng, 1)
  Ng(j,j) = 0;
end

geom.Neigs = Ng;

end


%% TEST ROUTINE
function test_geom_refneigs()

geom = geom_2dtri_longhex(5, 5, 20, 2, 'dir');
geom = geom_create_vacancies(geom, 1:11);
geom = geom_analyze(geom);
geom = geom_refneigs(geom, 2);

for n = 1:5:geom.nX
  In = find(geom.Neigs(:, n));
  triplot(geom.T', geom.X(1,:)', geom.X(2,:)'); view(2); hold on;
  plot(geom.X(1, n), geom.X(2,n), 'b*');
  plot(geom.X(1,In), geom.X(2, In), 'r*');
  hold off;
  pause;
end

end