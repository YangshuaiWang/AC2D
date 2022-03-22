
%% function [E, dE] = assemble_cb_toyeam(U, E, dE, geom, model)
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


function [E, dE] = assemble_cb_toyeam(U, E, dE, geom, model)

if nargin == 0
  test_assemble_cb_toyeam();
  return;
end

% extract model parameters for simpler code
a = model.a;
b = model.b;
c = model.c;
% rCutH = model.rCutH;  % >> not needed as it is encoded in neighbour
% matrix

% fixed model parameters
rho0 = 6 * exp(- b * 0.9);
rCB = model.rCB;

% read input
dDim = 2;
rDim = model.rDim;
U = reshape(U, rDim, geom.nX);
if nargout > 1
  dE = reshape(dE, rDim, geom.nX);
end

% loop over all elements
DW = zeros(2);
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
  r = Du * rCB;
  s = sqrt(r(1,:).^2 + r(2,:).^2);
  exp_a = exp(-a*(s-1)); exp_b = exp(-b*s);
  rho = sum( exp_b );
  W = sum( 0.5 * (exp_a.^2 - 2 * exp_a) ) ...
    + c * ( (rho - rho0)^2 + (rho-rho0)^4 );
  % NOTE: W is not rescale by det(A0). 
  % the volume is already corrected in geom.volT
  E = E + geom.volT(k) * W;
  
  % evaluate the CB forces
  if nargout > 1
%     [W, dW] = model.Wfun(model, Du);
%     % update total energy
%     E = E + geom.volT(k) * W;
%     % update gradient
%     dW = reshape(dW, rDim, dDim);
%     % dW taken w.r.t. Duhat instead of Du
%     dWhat = dW / J';
    dphi = -a * ( exp_a.^2 - exp_a ) ;
    dG = c * (2 * (rho-rho0) + 4 * (rho-rho0)^3);
    for i = 1:2
      for j = 1:2
        DW(i,j) = (  sum( (dphi./s) .* r(i,:) .* rCB(j,:) ) ...
          - b * dG * sum( (exp_b ./ s) .* r(i,:) .* rCB(j,:) )  );
      end
    end
    DW = DW / J';
    % write dWhat into dE
    for j = 1:2
      dE(:, t(j+1)) = dE(:, t(j+1)) + geom.volT(k) * DW(:,j);
      dE(:, t(1)) = dE(:,t(1)) - geom.volT(k) * DW(:,j);
    end
  end
end

% put into column vector (only required for testing!)
dE = dE(:);


end




%% test routine
% this routine trusts the implementation of E, but checks 
% the implementation of dE using finite differences.
function test_assemble_cb_toyeam()

% define model and geometry
a = 1; b = 3; c = 10; rCutH = 1;
model = model_toyeam_h(a, b, c, rCutH);
geom = geom_2dtri_mcrack(0, 4, 10, rCutH);
geom = bqc_prep_geom(geom, 2, 2, 'linearH');

disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);

% energy functional
E = 0;
dE = zeros(2, geom.nX);
fcnl = @(U_)(assemble_cb_toyeam(U_, E, dE, geom, model));
fcnl_old = @(U_)(assemble_cb(U_, E, dE, geom, model));

% base point
U = geom.X + 0.01 * rand(2, geom.nX);

[E, dE] = fcnl(U);
[E_old, dE_old] = fcnl_old(U);
abs(E - E_old)
norm(dE(:) - dE_old(:), inf)

% call finite difference test
addpath ./packages/popt
disp('--------------------------------------------');
disp('   Testing dE assembly in assemble_cb.m');
disp('--------------------------------------------');
test_derivatives(fcnl, U, 1);
disp('--------------------------------------------');

end



