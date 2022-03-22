
%% function nrm = W1p_norm(X, T, U, p)

function nrm = W1p_norm(X, T, U, p)

dDim = size(X, 1);
rDim = size(U, 1);

% loop over all elements
nrm = 0; 
fact_dDim = factorial(dDim);
for k = 1:size(T, 2)
  
  % compute Du
  t = T(:, k);
  J = zeros(dDim, dDim);
  Du = zeros(rDim, dDim);
  for j = 1:dDim
    J(:,j) = X(:,t(j+1)) - X(:,t(1));
    Du(:,j) = U(:,t(j+1)) - U(:,t(1));
  end
  if det(J) < 1e-5
    error('ERROR: small element in W1inf_norm');
  end
  % Du/Dx = Du/Dxhat * Dxhat/Dx = Du * J^{-1}
  Du = Du / J;
  
  % norm
  if p == inf
    nrm = max(nrm, max(abs(Du(:))));
  elseif 1 <= p
    nrm = nrm + det(J)/fact_dDim * sum(abs(Du(:)).^p);
  else
    error('ERROR: W1p_norm, p must be between 1 and inf');
  end
end

if p < inf
  nrm = nrm^(1/p);
end

end