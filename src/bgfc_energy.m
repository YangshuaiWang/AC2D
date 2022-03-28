%% function [E, dE] = bgfc_energy(U, geom, model)
%
% assembles simple energy and gradient for Blended Ghost Force Correction methods. Requires that
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

function [E, dE, hE, Ereal] = bgfc_energy(U, geom, model, useC, useP2)

if nargin == 0
  test_bgfc_energy();
  return;
end

if nargin < 4
  useC = false;
  useP2 = true;
end

% Prepare geom for assembly: 
% we assume that ***_prep_geom has already been called and put all
% important information into geom > there is nothing to do here.

% first convert the model
switch nargout
    case 1
        E = 0;
        % 1. assembly of atomistic component
        % 2. assembly of continuum component
        if useP2 == 0
%             E = assemble_a(U, E, [], geom, model, false, true);
%             E = assemble_cb(U, E, [], geom, model, true);
%             E = assemble_a(U, E, [], geom, model, false, false);
%             E = assemble_cb(U, E, [], geom, model, false);
            E = assemble_a(U, E, [], geom, model);
            E = assemble_cb(U, E, [], geom, model);
        else
%             E = assemble_a(U, E, [], geom, model, false, useC);
%             E = assemble_cb_p2(U, E, [], geom, model, false, useP2);
            E = assemble_a(U, E, [], geom, model);
            E = assemble_cb_p2(U, E, [], geom, model);
        end
        dE0 = model.dE0;
        dE0 = reshape(dE0, size(geom.X));
        U = reshape(U, size(geom.X));
        cor = dE0.*(U - model.F0*geom.X);
        E = E - sum(cor(:));        
    case 2
        E = 0; dE = zeros(size(U));
        % 1. assembly of atomistic component
        
        % 2. assembly of continuum component
        if useP2 ==0
            [E, dE] = assemble_a(U, E, dE, geom, model);
            [E, dE] = assemble_cb(U, E, dE, geom, model);
        else
%             [E, dE] = assemble_a(U, E, dE, geom, model);
%             [E, dE] = assemble_cb_p2(U, E, dE, geom, model);
        end
        % turn into column vector (required only for testing!)
        dE0 = model.dE0;
        dE0 = reshape(dE0, size(geom.X));
        U = reshape(U, size(geom.X));
        cor = dE0.*(U - model.F0*geom.X);
        E = E - sum(cor(:));
        dE = dE(:) - model.dE0;
        %E = E - su
        
%     case 3 
%         E = 0; dE = zeros(size(U));
%         % 1. assembly of atomistic component
%         [E, dE] = assemble_a(U, E, dE, geom, model);
%         % 2. assembly of continuum component
%         [E, dE] = assemble_cb_p2(U, E, dE, geom, model);
%         % global hessian
% %         hE = hEa + hEc;
%         % turn into column vector (required only for testing!)
%         dE0 = model.dE0;
%         dE0 = reshape(dE0, size(geom.X));
%         U = reshape(U, size(geom.X));
%         cor = dE0.*(U - model.F0*geom.X);
%         E = E - sum(cor(:));
%         dE = dE(:) - model.dE0;
end

% % OLD MATLAB VERSION
% switch nargout
%   case 1
%     E = 0;
%     % 1. assembly of atomistic component
%     E = assemble_a(U, E, [], geom, model);
%     % 2. assembly of continuum component
%     E = assemble_cb(U, E, [], geom, model);
%   case 2
%     E = 0; dE = zeros(size(U));
%     % 1. assembly of atomistic component
%     [E, dE] = assemble_a(U, E, dE, geom, model);
%     % 2. assembly of continuum component
%     [E, dE] = assemble_cb(U, E, dE, geom, model);
%     
%     % turn into column vector (required only for testing!)
%     dE = dE(:);  
% end

end


%% TEST ROUTINE
function test_bgfc_energy()
% define model and geometry
rCutH = 1; 
model = model_toyeam_h(4.0, 3.0, 10, rCutH);
% model = model_toyeam_h(4.0, 3.0, 10, 6*exp(-0.9*3.0), rCutH);
% geom = geom_2dtri_hexagon(20, 8, 2);
geom = geom_2dtri_mcrack(0, 3, 9, rCutH, 1.2);  
geom = geom_analyze(geom, rCutH);
geom0 = bgfc_prep_geom(geom, rCutH, 2, 'min-hessian', false);
geom = bgfc_prep_geom(geom, rCutH, 2, 'min-hessian', true);

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);
% base point
F0 = [1, 0; 0, 1]; model.F0 = F0;
% U = (F0+[0.01 -0.01; 0 0.01])*geom.X + 0.0 * rand(model.rDim, geom.nX);
U = (F0+0*[0.01 -0.01; 0 0.01])*geom.X;
U0 = (F0+0*[0.01 -0.01; 0 0.01])*geom0.X; 

useC = false; useP2 = false;
[E0, dE0] = bqce_energy(model.F0*geom0.X, geom0, model);
dE0 = reshape(dE0, size(geom0.X));
for n = 1: geom0.nX
    if geom0.di(n) > 2*rCutH + 1
        dE0(:, n) = 0;
    end
end
dE0(:,geom0.iBdry) = 0;
model0 = model; model0.E0 = E0; model0.dE0 = dE0(:);
temp = zeros(2, geom.nX); temp(:, 1:geom0.nX) = dE0; temp = temp(:);
model.E0 = E0; model.dE0 = temp;

% compare against Matlab code
[Ec, dEc] = bgfc_energy(U, geom, model); 

% ghost-force test
dEc = reshape(dEc, model.rDim, geom.nX);
dEc(:,geom.iBdry) = 0;
nrmG = norm(dEc(:), inf);
disp(['Ghost force test: ||dE(yA)||_inf = ', num2str(nrmG)]);
g = sum(abs(dEc), 1);
figure
trisurf(geom.T(1:3,:)', geom.X(1,:)', geom.X(2,:)', g');


format long
disp('Ghost Force test');
% disp(['   |err(E)| = ', num2str(abs(Ec - Em))]);
% disp(['   |err(dEm)| = ', num2str(norm(dEm, inf))]);
disp(['   |err(dEc)| = ', num2str(norm(dEc, inf))]);


% energy functional
useC = false; useP2 = true;
fcnl = @(U_)(bgfc_energy(U_, geom, model, useC, useP2));
useC = false; useP2 = false;
fcnl0 = @(U_)(bgfc_energy(U_, geom0, model0, useC, useP2));

% base point
rng('default')
U = (F0+[0.01 -0.01; 0 0.01])*geom.X + 0.02 * rand(model.rDim, geom.nX);
U0 = (F0+[0.01 -0.01; 0 0.01])*geom0.X + 0.02 * rand(model.rDim, geom0.nX);

% call finite difference test
addpath ./popt
disp('--------------------------------------------');
disp('   Testing dE assembly in bgfc_energy.m');
disp('--------------------------------------------');
test_derivatives(fcnl0, U0(:), 1);
disp('--------------------------------------------');


end

