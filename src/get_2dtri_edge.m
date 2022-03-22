
% get the edge dislocation for triangular lattice
% from JuLIPMaterials

function Xnew = get_2dtri_edge(X)


if nargin == 0
  Xnew = test_get_2dtri_edge();
  return;
end

xicorr = true;

Xnew = edge_predictor(X, xicorr);

end


%%%
%%% standard isotropic CLE edge dislocation solution
%%%
function Xnew = ulin_edge_isotropic(X)
    
    b = 1.0;
    nu = 0.25;
    x = X(1, :);
    y = X(2, :);
    r2 = x.^2 + y.^2;
    ux = b/(2*pi) * ( angle(x + sqrt(-1)*y) + (x .* y) ./ (2*(1-nu) * r2) );
    uy = -b/(2*pi) * ( (1-2*nu)/(4*(1-nu)) * log(r2) + - 2 * y.^2 ./ (4*(1-nu) * r2) );
    Xnew = [ux; uy];
end
% 

function xi1_ = xi1(x, y)
    b = 1.0;
    xi1_ = x - b * angle(x + sqrt(-1) * y) / (2*pi);
end

function dxi1_ = dxi1(xx, yy)
    b = 1.0;
    dxi1_ = 1 + b * yy / (xx^2 + yy^2) / (2*pi);
end

%%%
%%% lattice corrector to CLE edge solution; cf EOS paper
%%%
function xixixi = xi_solver(Y)
    Tol = 1e-10;
    maxnit = 10;
    yyy = Y(2);
    xxx = yyy;
    for n = 1:maxnit
        f = xi1(xxx, yyy) - Y(1);
        if abs(f) <= Tol
%             disp('newton solver done')
            break
        end
        xxx = xxx - f / dxi1(xxx, yyy);
    end
    if abs(xi1(xxx, yyy) - Y(1)) > Tol
        disp('newton solver did not converge at Y = $Y; returning input')
    end
    xixixi = [xxx, yyy];
end


%%%
%%% EOSShapeev edge dislocation solution
%%%
function XX = ulin_edge_eos(X)
    Xcr = zeros(2, size(X, 2));
    for n = 1:size(X,2)
        Xcr(:, n) = xi_solver(X(1:2,n));
    end
    XX = ulin_edge_isotropic(X);
end

% 
function Xnew = edge_predictor(X, xicorr)
   if xicorr
      Xnew = X(:,2:end) + ulin_edge_eos(X(:,2:end));
   else
      Xnew = X(:,2:end) + ulin_edge_isotropic(X(:,2:end));
   end
end


function Xnew = test_get_2dtri_edge()
close all
K = 5;


% lattice directions, auxiliary operators
a1 = [1;0];
Q6 = [cos(pi/3), -sin(pi/3); sin(pi/3), cos(pi/3)];
a2 = Q6 * a1;
a3 = Q6 * a2;
% aa = [a1, a2, a3, -a1, -a2, -a3, a1, a2];

%% atomistic core
% create a hexagon with 2*K+1 atoms across
X = zeros(2, 1 + 6*sum(1:K));
X(:,1) = [0;0]; ind = 1;
for j = 1:K
  x = Q6' * (j * a1 * ones(1, j) + a3 * (0:(j-1)));
  for k = 1:6
    x = Q6 * x;
    X(:, ind+1:ind+j) = x;
    ind = ind + j;
  end
end

Xnew = edge_predictor(X, true);
% Xnew = edge_predictor(X, false);

scatter(Xnew(1,:), Xnew(2,:))

end

