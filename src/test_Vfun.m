

function test_Vfun(model)

r = model.rCB + 0.01 * rand(size(model.rCB));
[V, dV] = model.Vfun(model, r);

for p = 2:10
  h = 10^(-p);
  dVh = zeros(size(dV));
  for j = 1:numel(r)
    rh = r; rh(j) = r(j) + h;
    Vh = model.Vfun(model, rh);
    dVh(j) = (Vh - V) / h;
  end
  
  disp(['p = ', num2str(p), '; err = ', num2str(norm(dV-dVh, inf))]);
%   disp(r(:,8:end)); disp(dV(:,8:end) - dVh(:,8:end));
%   pause;
end

end