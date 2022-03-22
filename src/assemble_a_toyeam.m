
%% function [E, dE] = assemble_a_toyeam(U, E, dE, geom, model, full_a)
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

function [E, dE] = assemble_a_toyeam(U, E, dE, geom, model)

if nargin == 0
  test_assemble_a_toyeam();
  return;
end

% this code ASSUMES that model.id = '2dtri_toyeam_h';
% it doesn't check whether it is.

% extract model parameters for simpler code
a = model.a;
b = model.b;
c = model.c;
% rCutH = model.rCutH;  % >> not needed as it is encoded in neighbour
% matrix

% fixed model parameters
rho0 = 6 * exp(- b * 0.9);

% reshape fields
U = reshape(U, model.rDim, geom.nX);
if nargout > 1
  dE = reshape(dE, model.rDim, geom.nX);
end

% the neighbour list in this version is static
Neigs = geom.Neigs;

%% loop over atom sites
for n = 1:geom.nX
  % if it isn't a proper atomistic site, skip it
  if geom.volX(n) <= 0
    continue;
  end
  
  % compute interaction neighbourhood 
  IN = find(Neigs(:, n));
  r = ( U(:, IN) - U(:, n) * ones(1,length(IN)) );
  s = sqrt(r(1,:).^2 + r(2,:).^2);
  
  % Assemble the energy
  ea = exp(-a*(s-1)); eb = exp(-b*s);
  phi = 0.5 * (ea.^2 - 2 * ea);
  E = E + geom.volX(n) * sum(phi);
  if c ~= 0
    rho = sum( eb );
    E = E + geom.volX(n) * c * ( (rho - rho0)^2 + (rho-rho0)^4 );
  end
  
  % Assemble the gradient
  if nargout > 1
    % pair interaction part
    dphi = - a * ( ea.^2 - ea );
    dV = ([1;1] * (dphi./s)) .* r;
    % eam part
    drho = - b * eb;
    dF = c * (2 * (rho - rho0) + 4 * (rho - rho0)^3);
    dV = dV + dF * ([1;1] * (drho./s)) .* r;
    % add to the global gradient
    dE(:, IN) = dE(:, IN) + geom.volX(n) * dV;
    dE(:, n) = dE(:, n) - geom.volX(n) * sum(dV, 2);
  end    
%       for i1 = 1:model.rDim
%         for i2 = 1:length(IN)
%           dE(i1, dof(IN(i2))) = dE(i1, dof(IN(i2))) + geom.volX(n) * dV(i1, i2);
%         end
%         dE(i1, dof(n)) = dE(i1, dof(n)) - geom.volX(n) * sum(dV(i1,:));
%       end
end

dE = dE(:);

end


%% TEST ROUTINE ASSEMBLY
function test_assemble_a_toyeam()
% define model and geometry
model = model_toyeam_h(4.0, 3.0, 5, 2);
geom = geom_2dtri_mcrack(3, 5, 5, 2);
geom.volX = ones(1,geom.nX);
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

