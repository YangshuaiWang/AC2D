%% function [E, dE] = bqc_energy(U, geom, model)
%
% assembles energy and gradient for blended QC methods. Requires that
% geom contains the fields volT, volX
%
% Input
%       U : (d x nX) nodal vector, deformation
%    geom : valid geometry structure, containing volT, volX
%   model : valid model structure (see model_toyeam, model_morse)
%
% Output 
%       E : energy
%      dE : gradient with respect to nodal values, without applying
%           the boundary conditions
%

function [E, dE] = bqc_energy(U, geom, model)

if nargin == 0
  test_bqc_energy();
  return;
end

% Prepare geom for assembly: 
% we assume that ***_prep_geom has already been called and put all
% important information into geom > there is nothing to do here.

switch nargout
  case 1
    E = 0;
    % 1. assembly of atomistic component
    E = assemble_a(U, E, [], geom, model);
    % 2. assembly of continuum component
    E = assemble_cb(U, E, [], geom, model);
  case 2
    E = 0; dE = zeros(size(U));
    % 1. assembly of atomistic component
    [E, dE] = assemble_a(U, E, dE, geom, model);
    % 2. assembly of continuum component
    [E, dE] = assemble_cb(U, E, dE, geom, model);
    
    % turn into column vector (required only for testing!)
    dE = dE(:);  
end

end


%% TEST ROUTINE
function test_bqc_energy()
% define model and geometry
model = model_morse(4.0, 1.0, 2);
geom = geom_2dtri_hexagon(10, 5);
geom = geom_analyze(geom);
geom = bqc_prep_geom(geom, 2, 2);

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);

% energy functional
fcnl = @(U_)(bqc_energy(U_, geom, model));

% base point
U = geom.X + 0.01 * rand(model.rDim, geom.nX);

% call finite difference test
addpath ./popt
disp('--------------------------------------------');
disp('   Testing dE assembly in assemble_a.m');
disp('--------------------------------------------');
test_derivatives(fcnl, U, 1);
disp('--------------------------------------------');

end

