
%% function [E, dE] = gqc23_energy(U, geom, model)
%
% assembles energy and gradient for for the GQC method. This method only
% implements the 2/3-method is implemented, which applies to
% nearest-neighbour multibody interactions in the triangular lattice.
%
% Input
%       U : (d x nX) nodal vector, deformation
%    geom : valid geometry structure runthrough geom_analyze
%   model : valid model structure, containing Vfun
%
% Output 
%       E : energy
%      dE : gradient with respect to nodal values, without applying
%           the boundary conditions
%

function [E, dE] = gqc23_energy(U, geom, model)

if nargin == 0
  test_gqc23_energy();
  return;
end

% check input: make sure it is a NN model
if model.rCutH > 1
  error('ERROR: gqc23_energy requires model.rCutH == 1');
end

% call different assembly components
switch nargout
  case 1
    E = 0;
    % 1. assembly of atomistic component
    E = assemble_a_toyeam(U, E, [], geom, model);
    % 2. assembly of continuum component
    E = assemble_cb_toyeam(U, E, [], geom, model);
    % 3. assembly of interface component
    E = assemble_gqc23(U, E, [], geom, model);
  case 2
    E = 0; dE = zeros(size(U));
    % 1. assembly of atomistic component
    [E, dE] = assemble_a_toyeam(U, E, dE, geom, model);
    % 2. assembly of continuum component
    [E, dE] = assemble_cb_toyeam(U, E, dE, geom, model);
    % 3. assembly of interface component
    [E, dE] = assemble_gqc23(U, E, dE, geom, model);

    % reshape gradient (needed mostly for testing)
    dE = dE(:);
end

end


%% main assembly routine
function [E, dE] = assemble_gqc23(U, E, dE, geom, model)

% resize dE
if nargout > 1, dE = reshape(dE, model.rDim, geom.nX); end
% preallocate some arrays
C = zeros(1,6);
r = zeros(2, 6);
% define modified modulus function
mod6 = @(i_)(mod(i_-1, 6)+1);  

% main loop for assembly of interface nodes
for n = 1:geom.nX
  % ignore if not an interface node
  if geom.volX(n) ~= -1
    continue;
  end
  
  % compute the 6 C parameters
  IN = sorted_neigs_2dtri(n, geom);
  for j = 1:6
    % if IN(j) is an interface or atomistic node then C(j) = 1.
    if abs(geom.volX(IN(j))) == 1
      C(j) = 1;
    % otherwise C = 2/3.
    else
      C(j) = 2/3;
    end
  end
  
  % compute the physical atom positions (relative to central atom)
  g = ( U(:, IN) - U(:, n) * ones(1,length(IN)) );
  % compute reconstructed positions (relative to central atom)
  for j = 1:6
    r(:, j) = C(j) * g(:, j) ...
      + (1-C(j)) * (g(:, mod6(j-1)) + g(:, mod6(j+1)));
  end

  % compute energy and gradient
  switch nargout
    case 1
      % compute V and add it to E
      V = model.Vfun(model, r);
      E = E + V;
      
    case 2
      % compute V and dV
      [V, dV] = model.Vfun(model, r);
      % add V to E
      E = E + V;
      % update dE
      for j = 1:6
        t = C(j) * dV(:,j) ...
          + (1-C(mod6(j-1))) * dV(:, mod6(j-1)) ...
          + (1-C(mod6(j+1))) * dV(:, mod6(j+1));
        dE(:, IN(j)) = dE(:, IN(j)) + t;
        dE(:, n) = dE(:, n) - t;
      end
  end
end

end


%% TEST ROUTINE
function test_gqc23_energy()
% define model and geometry
model = model_toyeam_h(4.0, 3.0, 10, 1); 
geom = geom_2dtri_mcrack(0, 5, 20, 1);
geomb = geom;
geom = gqc23_prep_geom(geom);
geomb = bqc_prep_geom(geomb, 1, 1, 'linearH');

geom.plot(geom);

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);

% make the ghost-force-test:
U = (eye(2) + 0.02 * rand(2)) * geom.X;
model.bc = 'dir';
[E, dE] = gqc23_energy(U, geom, model); 
% energy-test
Eb = bqc_energy(U, geomb, model);
disp(['Energy-Test: |Eq - Eb| / |Eb| = ', num2str(abs(E - Eb) / abs(Eb))]);
% ghost-force test
dE = reshape(dE, model.rDim, geom.nX); 

dE(:, 1:18) = 0;
dE(:, geom.iBdry) = 0;
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', abs(dE(1,:)));

nrmG = norm(dE(:), inf);
disp(['Ghost force test: ||dE(yA)||_inf = ', num2str(nrmG)]);

% energy functional
model.bc = '';
fcnl = @(U_)(gqc23_energy(U_, geom, model));

% base point
U = geom.X + 0.01 * rand(model.rDim, geom.nX);

% call finite difference test
addpath ./packages/popt
disp('--------------------------------------------');
disp('   Testing dE assembly in gqc23_energy.m');
disp('--------------------------------------------');
test_derivatives(fcnl, U, 1);
disp('--------------------------------------------');

end




