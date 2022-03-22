
%% function [E, dE] = assemble_cb(U, E, dE, geom, model)
%
% assembles the pure continuum components of a QC energy and adds to E, dE
%   E = E + \sum_{k} volT(k) W(Dy_k),   
%  dE = dE + associated gradient
%
% Input
%    U (nodal values vector), E (energy), dE (gradient),
%    geom (geometry structure), model (model structure)
% Output
%    E, dE are returned with updated values
%


% COMMENT: this routine does not know whether U is the displacement or
%          the deformation. This knowledge is only required by the
%          model structure!


function [E, dE] = assemble_cb(U, E, dE, geom, model)

if nargin == 0
  test_assemble_cb();
  return;
end

% read input
dDim = geom.dDim;
rDim = model.rDim;
U = reshape(U, rDim, geom.nX);
if nargout > 1
  dE = reshape(dE, rDim, geom.nX);
end

% loop over all elements
for k = 1:geom.nT
  % skip elements with zero effective volume
  if geom.volT(k) == 0
    continue;
  end
  
  % compute Du
  t = geom.T(:, k);
  J = zeros(dDim, dDim);
  Du = zeros(rDim, dDim);
  for j = 1:dDim
    J(:,j) = geom.X(:,t(j+1)) - geom.X(:,t(1));
    Du(:,j) = U(:,t(j+1)) - U(:,t(1));
  end
  if det(J) < 0.01
    error('ERROR: small element in assemble_cb');
  end
  % Du/Dx = Du/Dxhat * Dxhat/Dx = Du * J^{-1}
  Du = Du / J;
  
  % evaluate the Cauchy--Born energy
  switch nargout
    case 1
      W = model.Wfun(model, Du);
      E = E + geom.volT(k) * W;
      
    case 2
      [W, dW] = model.Wfun(model, Du);
      % update total energy
      E = E + geom.volT(k) * W;
      % update gradient
      dW = reshape(dW, rDim, dDim);
      % dW taken w.r.t. Duhat instead of Du
      dWhat = dW / J';
      % write dWhat into dE
      for j = 1:dDim
        dE(:, t(j+1)) = dE(:, t(j+1)) + geom.volT(k) * dWhat(:,j);
        dE(:, t(1)) = dE(:,t(1)) - geom.volT(k) * dWhat(:,j);
      end
      
    otherwise
      error('ERROR: assemble_cb must have 1 or 2 output arguments');
  end
end

% put into column vector (only required for testing!)
dE = dE(:);


end




%% test routine
% this routine trusts the implementation of E, but checks 
% the implementation of dE using finite differences.
function test_assemble_cb()

% define model and geometry
model = model_morse(4.0, 1.0, 2.2);
geom = geom_2dtri_hexagon(10, 6);
geom = geom_analyze(geom);
geom = bqc_prep_geom(geom, 2, 2);

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);

% energy functional
E = 0;
dE = zeros(2, geom.nX);
fcnl = @(U_)(assemble_cb(U_, E, dE, geom, model));

% base point
U = geom.X + 0.01 * rand(2, geom.nX);

% call finite difference test
addpath ./popt
disp('--------------------------------------------');
disp('   Testing dE assembly in assemble_cb.m');
disp('--------------------------------------------');
test_derivatives(fcnl, U, 1);
disp('--------------------------------------------');

end



