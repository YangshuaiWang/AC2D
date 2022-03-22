
function F = get_cb_groundstate(model)

if nargin == 0
  test_get_cb_groundstate();
  return;
end

F = eye(model.dDim);
[W, dW] = model.Wfun(model, F);
nit = 0;
a = 1;
theta = 0.01;
while ( (norm(dW(:)) > 1e-6) && (nit < 100) )
  % disp(norm(dW(:)))
  W1 = W; nrmdW2 = norm(dW(:))^2; a = a * 4;
  while ( (W1 > W - a * theta * nrmdW2) && (a > 1e-12) )
    a = a / 4;
    W1 = model.Wfun(model, F - a * dW);
  end
  a = a * 1.1;
  F = F - a * dW;
  [W, dW] = model.Wfun(model, F);
end


end

% function [W, dW] = Wwrap(F, model)
% [W, dW] = model.Wfun(reshape(F, rDim, dDim));
% dW = dW(:);
% end


function test_get_cb_groundstate()

model = model_toyeam(4, 3, 10);
F = get_cb_groundstate(model);

F

end