
% function geom = gqc23_prep_geom(geom)
%
% prepares geom so that it can be used to run gqc32_energy
% Note: this works only for nearest-neighbour models!
%
% Input
%     geom : valid geometry structure, triangular lattice.
%            If geom.volX exists, then with volX = 1 on nodes that
%            should be treated atomistically (including the interface!), 
%            and volX = 0 otherwise. If geom.volX does not exist, then
%            it is constructed to be maximal.
%
% Output
%   geom : * Modified volX: 1 in atomistic nodes, -1 for interface nodes,
%            0 elsewhere  (note that assemble_a will ignore nodes with
%            volX <= 0)
%          * Effective volumes of elements geom.volT
%

function geom = gqc23_prep_geom(geom)

if nargin == 0
  test_bqc23_prep_geom();
  return;
end

%% check validity of input and create preliminary volX
% if volX does not exist, construct default volX
if ~isfield(geom, 'volX')
  geom.volX = zeros(1, geom.nX);
  geom.volX(geom.di > 0) = 1;

% if it does exist already, make sure it is valid, i.e., all nodes
% inside the continuum region must be zero
else
  % first check that it is only 0s and 1s.
  if geom.volX ~= mod(geom.volX, 2)
    error('ERROR: gqc23_prep_geom; volX may contain only 0 and 1s.');
  end
  % now check that C-nodes are 0
  if sum(geom.volX(geom.di <= 0)) > 0
    error('ERROR: gqc23_prep_geom; volX must be 0 in continuum region.');
  end
end

%% create volT, and finalize volX
geom.volT = zeros(1, geom.nT);
facd = factorial(geom.dDim);
detA = det(geom.A);
for k = 1:geom.nT
  % calculate physical volume of T
  p = geom.X(:, geom.T(:, k));
  J = zeros(geom.dDim, geom.dDim);
  for j = 1:geom.dDim
    J(:, j) = p(:, j+1) - p(:, 1);
  end
  volT = abs(det(J)) / facd;
  
  % in Z^d each Voronoi cell has volume = 1, it has volume = det(A)
  % in A Z^d. But in the QC energy, we want to keep the effective
  % volume associated with each lattice site at 1, hence we need to
  % rescale the element volumes:
  volT = volT / detA;
  
  % add up the effective volumes of the vertices of T
  % > get number of vertices belonging to the atomistic region
  %  * WARNING: this part clearly works only for the 2D triangular lattice
  %  * note also that we have guaranteed above that if nA > 0 then T is a 
  %    nearest-neighbour triangle, hence the following works
  %  * the abs is because we set interface node volumes to -1
  nA = sum(abs(geom.volX(geom.T(:,k))));
  volT = volT * (1 - nA/3);
  
  % store volume of element
  geom.volT(k) = volT;

  % store interface node information:
  % if nA \in (0, 3) then at least one of the nodes must be atomistic and 
  % at one of them must be continuum and they are all nearest neighrbous
  % hence all atomistic nodes are actually interface nodes > set volumes
  % to -1.
  if ((nA > 0) && (nA < 3))
    for j = 1:3
      jX = geom.T(j,k);
      if geom.volX(jX) == 1
        geom.volX(jX) = -1;
      end
    end
  end
end

end


%% test routine
function test_bqc23_prep_geom()

geom = geom_2dtri_hexagon_2(40, 5, 3);
geom.fulla = false;
geom = geom_analyze(geom);
% to have one refined layer of pure CB, uncomment the next two lines.
% geom.volX = zeros(1, geom.nX);
% geom.volX(geom.di > 2) = 1;
geom = gqc23_prep_geom(geom);

geom.plot(geom);
for n = 1:geom.nX
  text(geom.X(1,n)+0.1, geom.X(2,n), num2str(geom.volX(n)));
end
for k = 1:geom.nT
  % get barycentre of element k
  x = 0.333 * sum(geom.X(:, geom.T(:, k)), 2);
  text(x(1), x(2), num2str(geom.volT(k), '%4.2f'));
end

end

