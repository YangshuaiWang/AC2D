
%% function [W, dW] = Wcb(model, F)
%
% Computation of the Cauchy--Born estored nergy density.
%
% Input: 
%        F : dxd matrix, or in 1 x d^2 vector in the format 
%            (e.g. in 2D) [F11,F21,F12,F22]
%            F is the deformation gradient
%                   (from the domain specified through the model)
%    model : model structure (e.g. taken from model_morse())
%
% Output:  W : scalar Wcb(F) (1x1)
%         dW : d x d matrix DWcb(F)
%

% TODO: write test routine (but it seems to work fine)

function [W, dW] = Wcb(model, F)
% read input
d = round(sqrt(numel(F)));
% reshape as matrix
F = reshape(F, d, d);

% Comments:
%  * we assume that the interaction neighbourhood is already given
%     in the model; this is called model.rCB, and is a subset of Z^d.
%  * the volume associated with a lattice site is just the cube,
%     i.e. volume = 1, hence W(F) = V( (F rCB_j)_j )

% Compute energy and gradient
if nargout == 1
  W = model.Vfun(model, F * model.rCB);
elseif nargout == 2
  [W, dV] = model.Vfun(model, F * model.rCB);
  % DW(F) = sum_j D_j V \otimes rCB_j   (easy to check)
  dW = dV * model.rCB';
else
  error('Wcb.m: Wcb must have 1 or 2 output arguments.');
end

end
