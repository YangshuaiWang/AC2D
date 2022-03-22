function [r, g] = get_grad_r(U, geom)

% read input
dDim = 2;
rDim = numel(U) / geom.nX;
U = reshape(U, rDim, geom.nX);

% preallocate r, g
r = zeros(1, geom.nX);
g = r;

% loop over all elements
for k = 1:geom.nT
  
  % compute Du
  t = geom.T(:, k);
  J = zeros(dDim, dDim);
  Du = zeros(rDim, dDim);
  for j = 1:dDim
    J(:,j) = geom.X(:,t(j+1)) - geom.X(:,t(1));
    Du(:,j) = U(:,t(j+1)) - U(:,t(1));
  end
  % Du/Dx = Du/Dxhat * Dxhat/Dx = Du * J^{-1}
  Du = Du / J;
  
  % evaluate r, g
  r(k) = norm(0.33 * sum(geom.X(:, t)));
  g(k) = norm(Du(:));
  
end

[r, I] = sort(r);
g = g(I);

end