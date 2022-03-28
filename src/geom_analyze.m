
%% function geom = geom_analyze(geom, rCutH)
%
% adds the following fields to the structure geom:
%  * nX, nT : number of vertices and number of elements)
%  * di     : array, di(n) = hopping distance of X(n) to the boundary of
%                    the refined region; negative in the continuum region
%  * NN : nearest neighbour matrix in delaunay triangulation 
%         e.g., use find(NN(:,n)) to get NNs of X(n) 
%  * geom_analyze = true
%  * E : edge list (useful for P2-fem, hessian assembly, etc)
%  * aTE, aET : adjacency relations
%  * iBdry : list of boundary nodes
%
function geom = geom_analyze(geom, rCutH)

% if nargin == 0 
%   test_geom_analyze(); 
%   return;
% end
if nargin == 0 
      test_geom_analyze(); 
      return;
elseif nargin == 1
      rCutH = 1;
end

% TERMINOLOGY: atomistic region = atomistically refined region
%                               = vertices where onL is true

% read data
T = geom.T';

%% 0. some simple things
geom.nX = size(geom.X, 2);
geom.nT = size(T, 1);


%% 1. create the nearest-neighbour relation as a sparse matrix
% TODO: implement this for d = 1, d = 3
switch geom.dDim
  case 2
    NN = sparse( [T(:,1); T(:,1); T(:,2); T(:,2); T(:,3); T(:,3)], ...
                 [T(:,2); T(:,3); T(:,1); T(:,3); T(:,1); T(:,2)], ...
                 0.1*ones(6*size(T,1), 1) );
    NN = ceil(NN);
  otherwise
    error('ERROR: geom_analyze is only implemented for d = 2 so far');
end
geom.NN = NN;

% compute full interaction neighbourhood
geom = geom_refneigs(geom, rCutH);


%% 2. determine the atomistic boundary
% case 1: geom_2dtri_hexagon, N == K
if ( geom.fulla )
  geom.di = geom.nX * ones(1, geom.nX);
  geom.di(geom.iBdry) = 0;
  I = geom.iBdry;

% case 2: anything else
else
  geom.di = (-1) * ones(1,geom.nX);
  I = [];
  for n = 1:geom.nX
    % if vertex is atomistic set di = nX, otherwise leave as -1 and continue
    if ~geom.onL(n)
      continue;
    end
    geom.di(n) = geom.nX;
    % get nearest neighbours
    iN = find(NN(:, n));
    % we now know that X(n) is in the atomistic region.
    % loop through neighbours and set di = 0 if one of the neighbours
    % are continuum; also store it in II
    for j = 1:length(iN)
      if ~geom.onL(iN(j))
        geom.di(n) = 0;
        I = [I, n];
        continue;
      end
    end
  end
end


%% 3. Calculate hopping distances
% for each atomistic vertex determine its hopping distance from the
% atomistic boundary
l = 0; 
while ~isempty(I)
  l = l+1; Inew = [];
  % loop through current "layer"
  % (note that the first layer is simply the atomistic boundary)
  for jI = 1:length(I)
    n = I(jI);
    % compute the neighours and give all atomistic neighbours with
    % di > l the value l and add to Inew
    iN = find(NN(:, n));
    for j = 1:length(iN)
      if ( geom.di(iN(j)) > l )
        geom.di(iN(j)) = l;
        Inew = [Inew, iN(j)];
      end
    end
  end
  I = Inew;
end


%% 4. do edge calculations
% also computes adjacency relations
geom = edges(geom);

% compute boundary (Dirichlet) nodes (boundary edges are iE)
if geom.fulla
  % for a full atomistic computation we need a buffer region
  % which is twice the interfaction radius
  geom.iBdry = find(geom.di <= 2*rCutH-1);
else
  % hence we get a complete list of boundary edges by taking the union
  % of all edge vertices
  iBdry = geom.E(:, geom.aET(2, :) == 0);
  % remove the "interior boundaries"
  iBdry( geom.di(iBdry) > 0 ) = [];
  % remove multiple elements
  geom.iBdry = unique(iBdry(:));
end


geom.geom_analyze = true;

end


%% function geom = edges(geom)
% Computes edges and adjacency relations for the mesh. This routine
% is only used for the initial mesh. Afterwards these relations are
% all automatically updated during refinement.
%
function geom = edges(geom)

% All edges with the interior edges duplicated
E = [geom.T([2,3], :), geom.T([1,3], :), geom.T([1,2], :)]; 
% sort each row so that the smaller index is first.
E = sort(E, 1); 
% throw away doubles and sort by rows(lexicographically)
[E, ~, J] = unique(E', 'rows'); E = E';
geom.E = E; geom.nE = size(geom.E, 2);

% find the adjacency relations between t and e. The order of the
% columns is such that (the local) edge j is opposite node j.
% The local basis has to be adjusted accordingly.
nT = geom.nT;
geom.aTE = [ J(1:nT), J((nT+1):(2*nT)), J((2*nT+1):(3*nT)) ]';

% Find the adjacency relation between E and T. The order
% defines also whether an element is considered "inside" the
% edge, or "outside" the edge. Boundary edges are always "inside".
% The mesh normal (DG) always points outside).
% *** VECTORIZE THIS *** (and also below in bisection())
geom.aET = zeros(2, size(E, 2));
for k = 1:nT
  t = geom.T(:, k); t = [t; t(1)];  %#ok
  for j = 1:3
    % case 1: left
    if geom.E(1, geom.aTE(j, k)) == t(j+1)
      geom.aET(1, geom.aTE(j, k)) = k;
    % case two: right
    else
      geom.aET(2, geom.aTE(j, k)) = k;
    end
  end
end

% Finally, fix the direction of the edges. For all boundary edges we need 
% to check whether they are oriented in such a way that the (unique)
% neighbouring element is to the left.
iE = find(geom.aET(1, :) == 0);
geom.E(:, iE) = geom.E([2,1], iE);
geom.aET(:, iE) = geom.aET([2,1], iE);


end




%% test-code for geom_analyze
%
function test_geom_analyze()
geom = geom_2dtri_hexagon(10, 5, 1, 'per');
geom = geom_create_vacancies(geom, [1,2, 30]);
geom = geom_analyze(geom);
geom.plot(geom);
for n = 1:geom.nX
  text(geom.X(1,n)+0.2, geom.X(2,n), num2str(geom.di(n)));
end
end

