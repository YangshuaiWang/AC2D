
%% function [E, dE] = bqc_energy(U, geom, model)
%
% assembles energy and gradient for the full atomistic model. Requires that
% geom.X, geom.T represents the full a reference lattice.
%
% Input
%       U : (d x nX) nodal vector, deformation
%    geom : valid geometry structure, containing volT, volX
%   model : valid model structure (see model_toyeam, model_morse)
%
% Output 
%       E : energy
%      dE : gradient with respect to nodal values, 
%           if bc == 'per', then the b.c. is already applied
%           and the superfluous components set to 0.
%

function [E, dE] = atm_energy(U, geom, model)

if nargin == 0
  test_atm_energy();
  return;
end

switch nargout
  case 1
    E = 0;
    E = assemble_a(U, E, [], geom, model, true);
  case 2
    E = 0; dE = zeros(size(U));
    [E, dE] = assemble_a(U, E, dE, geom, model, true);
    
    % turn into column vector (required only for testing!)
    dE = dE(:);  
end

end


%% TEST ROUTINE
function test_atm_energy()
% define model and geometry
N = 5;
model = model_toyeam_h(4.0, 2.0, 10, 1);
geom = geom_2dtri_hexagon(N, N, 1, 'dir');
geom = geom_create_vacancies(geom, 1);
geom = geom_analyze(geom);
geom.volX = ones(1, geom.nX);
model.F0 = [1, 0; 0, 1];

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);

% wrapper functional
  function [E, dE] = atm_wrapper(W, geom, model, iFree)
    Y = geom.X;
    Y(:, iFree) = reshape(W, 2, length(iFree));
    if nargout == 1
      E = atm_energy(Y, geom, model);
    else
      [E, dE] = atm_energy(Y, geom, model);
      dE = reshape(dE, 2, geom.nX);
      dE = dE(:, iFree);
      dE = dE(:);
    end
  end

% energy functional
% iFree = setdiff(1:geom.nX, geom.iPer(1,:));
if isfield(geom, 'iDir')
  iFree = setdiff(1:geom.nX, geom.iDir);
else
  iFree = setdiff(1:geom.nX, geom.iBdry);
end

fcnl = @(W_)(atm_wrapper(W_, geom, model, iFree));

% base point
U = geom.X + 0.01 * rand(model.rDim, geom.nX);
W0 = U(:, iFree);

% finite difference test
addpath ./packages/popt
disp('--------------------------------------------');
disp('   Testing dE assembly in assemble_a.m');
disp('--------------------------------------------');
test_derivatives(fcnl, W0(:), 1);
disp('--------------------------------------------');

end

