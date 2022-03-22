
function v = geom_volume(geom)

v = 0;
dDim = geom.dDim;

% loop over all elements
for k = 1:geom.nT
  % compute Du
  t = geom.T(:, k);
  J = zeros(dDim, dDim);
  for j = 1:dDim
    J(:,j) = geom.X(:,t(j+1)) - geom.X(:,t(1));
  end
  v = v + det(J)/2;
end

end
