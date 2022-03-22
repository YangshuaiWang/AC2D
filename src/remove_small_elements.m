
function T = remove_small_elements(T, X, d0)

% loop over all elements
kDEL = zeros(1, size(T,2));
for k = 1:size(T, 2)
  t = T(:, k);
  J = zeros(2, 2);
  diamT = 0;
  for j = 1:2
    J(:,j) = X(:,t(j+1)) - X(:,t(1));
    diamT = max(diamT, norm(J(:,j)));    
  end
  if det(J)/diamT < d0
    kDEL(k) = true;
  end
end
T(:, kDEL==true) = [];
end
