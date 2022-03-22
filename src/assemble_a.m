
%% function [E, dE] = assemble_a(U, E, dE, geom, model, full_a)
%
% assembles the pure atomistic components of a QC energy and adds to E, dE
%   E = E + \sum_{n} volX(n) E_n,   
%  dE = dE + associated gradient
%
% Input
%    U (nodal values), E (energy), dE (gradient),
%    geom (geometry structure), model (model structure)
%    full_a : optional flag, if true then it is assumed that a
%             full atomistic model is used, which requires a slightly
%             different treatment of the boundary
% Output
%    E, dE are returned with updated values
%


% COMMENT: this is an assembly with loops over all atoms, but it should
%          be easy to write a vectorized version. The problem though
%          is that a vectorized version might need an insane amount of
%          memory since all finite difference stencils need to be stored
%          at the same time. An alternative might be to rewrite
%          assemble_a in Fortran or C

% COMMENT: this routine does not know whether U is the displacement or
%          the deformation. This knowledge is only required by the
%          model structure!
%       >> actually, this is not quite true. In the current datastructure
%          the Vfun needs deformed positions, not displacements!
%       >> fix this in future versions?

function [E, dE] = assemble_a(U, E, dE, geom, model, full_a)

if nargin == 0
  test_assemble_a();
  % test_neighbour_matrix();
  % test_neighbour_matrix_fulla();
  return;
end

% if the model is the toy_eam model with static interaction
% then call a simpler assembly routine.
% if strcmp(model.id, 'tri2d_toyeam_h')
%   if nargout == 1
%     E = assemble_a_toyeam(U, E, dE, geom, model);
%   else
%     [E, dE] = assemble_a_toyeam(U, E, dE, geom, model);
%   end
%   return;
% end


if nargin < 6
  full_a = false;
end

% reshape fields
U = reshape(U, model.rDim, geom.nX);
if nargout > 1
  dE = reshape(dE, model.rDim, geom.nX);
end

% dof list: dof(i) = i  (changes only for periodic b.c.)
dof = 1:geom.nX;
% also write the periodicity information in a large array
if full_a && strcmp(model.bc, 'per')
  %  dof(i) = j   if the DOF i is associated with the DOF j
  dof(geom.iPer(1,:)) = geom.iPer(2,:);
  % arrange to skip all "multiple" nodes
  geom.volX = ones(1, geom.nX);
  geom.volX(geom.iPer(1,:)) = 0;
  
  aper = [ [0;0], model.F0 * geom.aper ];
  Lper = norm(geom.aper(:,1)) / 2;
end

%% Neighbour List
% in this version of the code, we assume that the neigbour list
% precomputed!
Neigs = geom.Neigs;
% Neigs = neighbour_matrix(U, geom, model, full_a);


%% loop over atom sites
for n = 1:geom.nX
  % if it isn't a proper atomistic site, skip it
  if geom.volX(n) <= 0
    continue;
  end
  
  % compute interaction neighbourhood 
  IN = find(Neigs(:, n));
  r = ( U(:, IN) - U(:, n) * ones(1,length(IN)) );
  
  
  % if we are in a fully atomistic problem and are close to the
  % domain boundary, then we should correct "r"
  % TODO: distinguish periodic and Dirichlet b.c.!
  if full_a && (geom.di(n) <= model.rCutH) && strcmp(model.bc, 'per')
    geom_per_reconstruct_c(r, aper, Lper);
  end


%   % DEBUGGING the interaction ranges!
%   geom.plot(geom); hold on;
%   plot(geom.X(1,n), geom.X(2,n), 'k*', 'Markersize', 20);
%   plot(geom.X(1, IN), geom.X(2,IN), 'k.', 'Markersize', 20); hold off;
%   pause;
  
  
  switch nargout
    case 1
      % compute V and add it to E
      V = model.Vfun(model, r);
      E = E + geom.volX(n) * V;
      
    case 2
      % compute V and dV
      [V, dV] = model.Vfun(model, r);
      % add V to E
      E = E + geom.volX(n) * V;
      % update dE
      dE(:, dof(IN)) = dE(:, dof(IN)) + geom.volX(n) * dV;
      dE(:, dof(n)) = dE(:, dof(n)) - geom.volX(n) * sum(dV, 2);
%       for i1 = 1:model.rDim
%         for i2 = 1:length(IN)
%           dE(i1, dof(IN(i2))) = dE(i1, dof(IN(i2))) + geom.volX(n) * dV(i1, i2);
%         end
%         dE(i1, dof(n)) = dE(i1, dof(n)) - geom.volX(n) * sum(dV(i1,:));
%       end
  end
end

dE = dE(:);

end





%% compute neighbour list of atomistic part of of deformed configuration
% using delaunay triangulation
function [Neigs, T, Z, Neigs0] = neighbour_matrix(U, geom, model, full_a)


% two cases: full_a = true / false
% if false, then we just triangulate the atomistic region and are done.
% for full_a we have to work a bit more

% WARNING: for general rDim, U doesn't give positions!!!
%          TODO: write a routine that reconstructs node positions from U

if ( isfield(geom, 'neig_useT') && (geom.neig_useT == true) )
  % use the reference triangulation
  T = geom.T';
  % prepare for reconstruction of boundary stencils
%   aper = [ [0;0], model.F0 * geom.aper ];
%   Lper = norm(geom.aper(:,1)) / 2;
elseif (~full_a) || (~ strcmp(model.bc, 'per'))
  % case 1: small central atomistic region
  Z = U(:, geom.onL == 1);
  DT = DelaunayTri(Z');
  T = DT.Triangulation;
  T = remove_small_elements(T', Z, 0.1)';
else
  % case 2: complete atomistic region
  % compute constrained delaunay triangulation:
  %   * nodes are interpolation between U and X positions
  %   * constraints are edges at the boundary
  N = geom.N;
  U = reshape(U, model.rDim, geom.nX);
  r = sqrt(sum(geom.X.^2, 1));
  lambda = (r / N - 3/4) / (1/2 - 3/4);
  lambda = [1;1] * min(max(lambda, 0), 1);
  Z = U .* lambda + geom.X .* (1-lambda);

  % delaunay
  C = geom.E(:, geom.aET(2,:) == 0);
  DT = DelaunayTri(Z', C');
  T = DT.Triangulation;
  T = remove_small_elements(T', Z, 0.1)';
  % prepare for reconstruction of boundary stencils
%   aper = [ [0;0], model.F0 * geom.aper ];
%   Lper = norm(geom.aper(:,1)) / 2;
end


% compute adjacancy matrix
% TODO: 2D only so far!
NN = sparse( [T(:,1); T(:,1); T(:,2); T(:,2); T(:,3); T(:,3)], ...
             [T(:,2); T(:,3); T(:,1); T(:,3); T(:,1); T(:,2)], ...
              ones(6*size(T,1), 1) );

% if the model is fully atomistic and the boundary conditions are periodic
% then adjust the NN matrix to connect periodically repeated images
if full_a && strcmp(geom.bc, 'per')
    for j = 1:size(geom.iPer, 2)
      NN(:, geom.iPer(2,j)) = NN(:,geom.iPer(2,j)) + NN(:,geom.iPer(1,j));  %#ok
    end
    % symmetrize
    NN = NN + NN';
    % remove interactions with "ghost atoms"
    NN(geom.iPer(1,:), :) = 0;
    NN(:, geom.iPer(1,:)) = 0;
end

% get 2nd, 3rd, ..., rCutH neighours as well
Neigs = NN;
for j = 2:model.rCutH
  Neigs = Neigs * NN;
end
% remove self-interaction  (slow!)
for n = 1:size(Neigs, 1)
  Neigs(n,n) = 0.1;   %#ok
end

% % remove all interactions outside the interaction range (very slow!)
% for n = 1:size(Neigs, 2)
%   I = find(Neigs(:, n));
%   dZ = Z(:, I) - Z(:, n) * ones(1,length(I));
%   if full_a && (geom.di(n) <= model.rCutH)
%     geom_per_reconstruct_c(dZ, aper, Lper);
%   end
%   dist2 = sum( dZ.^2, 1 );
%   Neigs(I(dist2 > model.rCut^2), n) = 0.1;
% end

Neigs = round(Neigs);

% write Neigs matrix into large matrix
if ~full_a
  I = find(geom.onL);
  Neigs0 = Neigs;
  Neigs = spalloc(geom.nX, geom.nX, nnz(Neigs0));
  Neigs(I, I) = Neigs0;
end


end






%% TEST ROUTINE ASSEMBLY
function test_assemble_a()
% define model and geometry
model = model_toyeam2(4.0, 3.0, 10);
geom = geom_2dtri_hexagon(10, 7, 2);
geom = geom_analyze(geom);
geom = bqc_prep_geom(geom, model.rCutH, 2, 'linearH');

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);

% energy functional
E = 0;
dE = zeros(2, geom.nX);
fcnl = @(U_)(assemble_a(U_, E, dE, geom, model));

% base point
U = geom.X + 0.01 * rand(2, geom.nX);

% call finite difference test
addpath ./popt
disp('--------------------------------------------');
disp('   Testing dE assembly in assemble_a.m');
disp('--------------------------------------------');
test_derivatives(fcnl, U, 1);
disp('--------------------------------------------');
end


%% TEST ROUTINE NEIGHBOUR MATRIX
function test_neighbour_matrix()

% define model and geometry
model = model_toyeam2(4.0, 3.0, 10);
geom = geom_2dtri_hexagon(20, 10, 2);
Iv = find( (abs(geom.X(2,:)) < 1/4) & (abs(geom.X(1,:)) < 3.5) );
geom = geom_create_vacancies(geom, Iv);
geom = geom_analyze(geom);
geom = bqc_prep_geom(geom, model.rCutH, 2, 'linearH');

U = geom.X + 0.01 * rand(2, geom.nX);
[Neigs, T, Z, Neigs0] = neighbour_matrix(U, geom, model, false);

% % correctness test 1
% n = 95;
% In = find(Neigs0(:, n));
% 
% trisurf(T, Z(1,:)', Z(2,:)', 0 * Z(1,:)'); view(2); hold on;
% plot(Z(1,n), Z(2,n), 'b*');
% plot(Z(1,In), Z(2, In), 'r*');
% hold off;

% correctness test 2
n = 300;
In = find(Neigs(:, n));

trisurf(T, Z(1,:)', Z(2,:)', 0 * Z(1,:)'); view(2); hold on;
plot(geom.X(1, n), geom.X(2,n), 'b*');
plot(geom.X(1,In), geom.X(2, In), 'r*');
hold off;

% % timing test
% profile on;
% tic
% for n = 1:10
%   Neigs = neighbour_matrix(U, geom, model, false);
% end
% toc
% profile off;
% profile report

end


%% TEST ROUTINE NEIGHBOUR MATRIX
function test_neighbour_matrix_fulla()

% define model and geometry
model = model_toyeam2(4.0, 3.0, 10);

K0 = 3;
geom = geom_2dtri_longhex(K0, 10, 10, 2, 'per');

% Iv = find( (abs(geom.X(2,:)) < 1/4) & (abs(geom.X(1,:)) < 3.5) );
% geom = geom_create_vacancies(geom, Iv);
geom = geom_create_vacancies(geom, 1:2*K0+1);

geom = geom_analyze(geom);
model.F0 = [1, 0.03; 0, 1.03];

U = 0.01 * rand(2, geom.nX);
U(geom.iPer(1,:)) = U(geom.iPer(2,:));
Y = model.F0 * geom.X + U;

% correctness test 2
[Neigs, T, Z] = neighbour_matrix(Y, geom, model, true);

for n = 1:10:geom.nX
  In = find(Neigs(:, n));
  triplot(T, Y(1,:)', Y(2,:)'); view(2); hold on;
  plot(Y(1, n), Y(2,n), 'b*');
  plot(Y(1,In), Y(2, In), 'r*');
  hold off;
  pause;
end

% % timing test
% profile on;
% tic
% for j = 1:10
%   E = assemble_a(U, 0, [], geom, model, true);
% end
% toc
% profile off;
% profile report

end
