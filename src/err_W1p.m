% [err1, err2, errinf] = err_W1p(U1, geom1, U2, geom2)
%
% approximately computes the W1p-norm between (U1, geom1) and (U2, geom2)
% as follows
%    err = max( || U1 - I1(U2) ||, || I2(U1) - U2 || )       (*)
% where Ij(Ui) is the geomj-interpolant of Ui, and ||.|| the W1p-norm
%
% INPUT
%    Ui, geomi : nodal values Ui on geometry geomi
%            p : Sobolev index for W^{1,p}-norm
% OUTPUT
%          err : approximation of W1p-distance; cf. (*)
%

function [err1, err2, errinf] = err_W1p(U, geom, U_ref, geom_ref)

if nargin == 0
  test_err_W1p();
  return;
end

U_int = eval12(U, geom, geom_ref);
err1 = W1p_norm(geom_ref.X, geom_ref.T, U_int - U_ref, 1);
err2 = W1p_norm(geom_ref.X, geom_ref.T, U_int - U_ref, 2);
errinf = W1p_norm(geom_ref.X, geom_ref.T, U_int - U_ref, inf);

end




%% function U2 = eval12(U1, geom1, geom2)
% evaluate (U1, geom1) in geom2 and return the new nodal values
function U2 = eval12(U1, geom1, geom2)

dDim = geom1.dDim;
rDim = size(U1, 1);
U2 = zeros(rDim, geom2.nX);

% get a delaunay triangulation of geom1.X,
% but we need to cheat a bit to make the domain numerically convex
% and to ensure that all points really lie inside (stretch by 1e-5!)
% but this is only used to create the triangulation, for the evaluation
% of U2 the original nodes are used
X = geom1.X; N = geom1.N;
r = sqrt(sum(X.^2, 1));
X(:, geom1.iBdry) = X(:, geom1.iBdry) .* ...
  ( ones(dDim, 1) * (1+1e-5 + 1e-1 * (N-r(geom1.iBdry))/N ) );
DT = DelaunayTri(X');
% find the elements in geom1 to which the vertices in geom2 belong
IT = pointLocation(DT, geom2.X');
% compute the local convex coordinates of these vertices
for i = 1:geom2.nX
  iT = IT(i);
  % if none was found, then set U to zero!
  if isnan(iT)
    U2(:, i) = 0;
    continue;
  end
  % compute Du1
  t = DT.Triangulation(iT, :);
  J = zeros(dDim, dDim);
  Du = zeros(rDim, dDim);
  for j = 1:dDim
    J(:,j) = geom1.X(:,t(j+1)) - geom1.X(:,t(1));
    Du(:,j) = U1(:,t(j+1)) - U1(:,t(1));
  end
  if det(J) < 0.01
    error(['Small element in err_W1p:eval12, J = ', num2str(det(J))]);
  end
  % Du/Dx = Du/Dxhat * Dxhat/Dx = Du * J^{-1}
  Du = Du / J;
  
  % evaluate U2
  U2(:, i) = U1(:, t(1)) + Du * (geom2.X(:, i) - geom1.X(:, t(1)));
end

end


%% test routine
function test_err_W1p()

% N = 30;
% 
% geom1 = geom_2dtri_hexagon_2(N, 5, 5);
% geom1 = geom_analyze(geom1);
% r1 = sqrt(sum(geom1.X.^2, 1));
% U1 = r1.^(-3); U1 = U1 + 1e-4 * rand(size(U1));
% % U1 = exp(geom1.X(1,:) .* geom1.X(2,:) / N^2);
% 
% geom2 = geom_2dtri_hexagon(N, 10, 1);
% geom2 = geom_analyze(geom2);
% r2 = sqrt(sum(geom2.X.^2, 1));
% U2 = r2.^(-3); U2 = U2 + 1e-4 * rand(size(U2));
% % U2 = exp(geom2.X(1,:) .* geom2.X(2,:) / N^2);
% 
% err = err_W1p(U1, geom1, U2, geom2, inf);
% disp(['Error = ', num2str(err)]);
% 
% U12 = eval12(U1, geom1, geom2);
% 
% figure(1);
% subplot(2,2,1);
% trisurf(geom1.T', geom1.X(1,:)', geom1.X(2,:)', U1);
% subplot(2,2,2);
% trisurf(geom2.T', geom2.X(1,:)', geom2.X(2,:)', U12);
% subplot(2,2,3);
% trisurf(geom2.T', geom2.X(1,:)', geom2.X(2,:)', U2 - U12);
% 

end