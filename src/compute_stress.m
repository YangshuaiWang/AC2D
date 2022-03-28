% function [Du, dW] = compute_stress(X, U, model)
% Program: for input element's coordinates X and deformation U, compute its
%       gradient Du and CB stress tensor according to energy func in model.
function [Du, dW] = compute_stress(X, U, model)
dDim = 2;
rDim = 2;
J = zeros(dDim, dDim);
Du = zeros(rDim, dDim);
% X = X(:,end:-1:1);
% U = U(:,end:-1:1);
for j = 1:dDim
    J(:,j) = X(:, j+1) - X(:, 1);
%     J(:,j) = geom.X(:,t(j+1)) - geom.X(:,t(1));
    Du(:,j) = U(:, j+1) - U(:, 1);
end
if det(J) < 0.01
    keyboard;
    error('ERROR: small element in get_model_error');
end
% Du/Dx = Du/Dxhat * Dxhat/Dx = Du * J^{-1}
Du = Du / J;
if nargout == 1;
    return;
end

% evaluate the Cauchy--Born stress tensor

[~, dW] = model.Wfun(model, Du);
% update gradient
dW = reshape(dW, rDim, dDim);
dW = 2*dW/sqrt(3);
end