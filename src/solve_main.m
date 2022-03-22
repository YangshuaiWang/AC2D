% function Y = solve_main(geom, model, F0)
%
% This routine takes on various tasks required to solve the nonlinear
% optimisation problem associated with QC methods:
%   * apply the boundary conditions. This is achieved by writing
%     Y = F0 * X + W, removing unnecessary dofs from the displacement W 
%     (e.g. dirichlet nodes, or duplicated nodes for periodic bcs)
%     and then wrapping the energy functionals in a wrapper that transforms
%     between U and reduced W.
%   * provide an interface to nonlinear and linear solvers
%   * provide various preconditioners (though not at the moment...)
%
% INPUT
%   geom, model : valid geometry and model structures
%            F0 : macroscopic strain s.t. Y = F0 * X + W. 
% OUTPUT
%             Y : nodal values of optimised deformation
%             E : energy
%
% TODO: * create more solver options
%       * move the geometry analysis, etc to this solver
%         and instead create and options structure for the a/c methods
%

function [Y, E] = solve_main(geom, model, F0, opts)

model.F0 = F0;

% check F0
if ( (size(F0, 1) ~= model.rDim) || (size(F0,2) ~= model.dDim) )
  error('ERROR: solve_main requires size(F0) = [rDim, dDim]');
end

% compute a list of free indices
switch model.bc
  case 'dir'
    if isfield(geom, 'iDir')
      iFree = setdiff(1:geom.nX, geom.iDir);
    else
      iFree = setdiff(1:geom.nX, geom.iBdry);
    end
  case 'per'
    iFree = setdiff(1:geom.nX, geom.iPer(1,:));
  case ''
    iFree = 1:geom.nX;
  otherwise
    error('ERROR: solve_main encountered an unknown bdry condition');
end

% create an energy functional
switch opts.actype
  case 'bqc'
    energy = @(Y_)(bqc_energy(Y_, geom, model));
  case 'gqc23'
    energy = @(Y_)(gqc23_energy(Y_, geom, model));
  case 'qce'
    % this assumes that the setup of geom is such that it
    % reduced to QCE
    energy = @(Y_)(bqc_energy(Y_, geom, model));
  case 'cb'
    energy = @(Y_)(cb_energy(Y_, geom, model));
  case 'scb'
    energy = @(Y_)(scb_energy(Y_, geom, model));
  case 'atm'
    energy = @(Y_)(atm_energy(Y_, geom, model));
  otherwise
    error('ERROR: solve_main encountered an unknown a/c method type');
end

% construct the preconditioner
switch opts.prec
  case 'laplace'
    P = prec_laplace(geom, model);
  case 'dyn_laplace'
    P = @(W_)( prec_dyn_laplace( ...
      get_dfm(W_, geom, model, F0, iFree), geom, model, [0.3, 1e-3] ) );
  case 'none'
    P = 1;
  otherwise 
end

% visualisation of optimisation
if isfield(opts, 'visualize') ...
    && opts.visualize && isempty(opts.popt.visualize)
  opts.popt.visualize = @(W_)(...
    visualize(get_dfm(W_, geom, model, F0, iFree), geom));
end

% initial condition
if isempty(opts.U0)
  W = zeros(model.rDim*length(iFree), 1);
else
  W = opts.U0 - F0 * geom.X;
  W = W(:, iFree);
  W = W(:);
end


% TODO: check that solver is popt, or different solver and call the
% right solver

% run nonlinear solver
addpath ./packages/popt
ffun = @(W_)(energy_wrapper(W_, geom, model, F0, iFree, energy));

% call popt
W = poptls(ffun, P, W, opts.popt);

% to improve the accuracy, run also fsolve
if opts.post_fsolve
  gfun = @(W_)(get_grad(W_, ffun));
  W = fsolve(gfun, W);
end

if nargout > 1, E = ffun(W); end
Y = get_dfm(W, geom, model, F0, iFree);

end

% extract gradient only
function dE = get_grad(u, ffun)
[~, dE] = ffun(u);
end

%% Y = get_dfm(W, geom, model. F0, iFree)
% recreate the deformation field from W and the boundary conditions
function Y = get_dfm(W, geom, model, F0, iFree)
U = zeros(model.rDim, geom.nX); 
U(:, iFree) = reshape(W, model.rDim, length(iFree));
if strcmp(model.bc, 'per')
  U(:, geom.iPer(1,:)) = U(:, geom.iPer(2,:));
end
Y = F0 * geom.X + U;
end


%% energy_wrapper
% applies the boundary conditions
function [E, dE] = energy_wrapper(W, geom, model, F0, iFree, energy)

% recreate the deformation field from W and the boundary conditions
Y = get_dfm(W, geom, model, F0, iFree);

% evaluate the energy and gradient
if nargout == 1
  E = energy(Y);
  
elseif nargout == 2
  [E, dE] = energy(Y);
  % apply the bc to the gradient
  dE = reshape(dE, model.rDim, geom.nX);
  switch model.bc
    case 'dir'
      dE = dE(:, iFree);
    case 'per'
      for j = 1:size(geom.iPer, 2)
        dE(:,geom.iPer(2,j)) = dE(:,geom.iPer(2,j)) + dE(:,geom.iPer(1,j));
      end
      dE = dE(:, iFree);
    % case ''
      % if there are no constraints, do nothing
  end
  % reshape gradient so that it can be used by the optimisation routine
  dE = dE(:);
end

end





function visualize(Y, geom)

hdl = trisurf(geom.T', Y(1,:)', Y(2,:)', 0 * Y(1,:)'); hold on;
set(hdl, 'Facecolor', 'white', 'Edgecolor', 'black');
plot(Y(1,:), Y(2,:), 'b.', 'Markersize', 10);
hold off;
view(2);
drawnow;

end


