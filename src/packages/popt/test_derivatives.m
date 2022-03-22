%
% TEST_DERIVATIVES
%
% [errG, errH] = test_derivatives(fun, u, ndiff)
%
% Uses finite differences to test the implementation of the functional
% fcnl which is of type [F, G, H] = fcnl(u).
% Returns the last error of the hessian
%
function [errG, errH] = test_derivatives(fun, u, ndiff)

if nargin < 3, ndiff = 2; end

if ndiff == 2
  [errG, errH] = test_diff_2(fun, u);
elseif ndiff == 1
  errG = test_diff_1(fun, u);
  errH = [];
else
  error('invalid ndiff input');
end

end


function [errG, errH] = test_diff_2(fun, u)

[F, G, H] = fun(u);
N = length(u);

disp('Testing the gradient');
for p = 2:12
    
    h = 10^(-p);
    Gh = zeros(N, 1);
    
    for i = 1:N
        uh = u; uh(i) = uh(i) + h;
        Fh = fun(uh);
        Gh(i) = (Fh - F) / h;
    end
    
    disp([ 'p = ', num2str(p), ', err = ', num2str(norm(G- Gh, inf)) ]);
end

errG = G - Gh;

disp(' ');
disp('Testing the Hessian');

for p = 2:10
    
    h = 10^(-p);
    Hh = zeros(N, N);
    
    for i = 1:N
        uh = u; uh(i) = uh(i) + h;
        [Fh, Gh] = fun(uh);
        Hh(:,i) = (Gh - G) / h;
    end
    
    error = full(H - Hh);
    error = max(abs(error(:)));
    disp([ 'p = ', num2str(p), ', err = ', num2str(error)]);
end


errH = H - Hh;

end



function errG = test_diff_1(fun, u)

[F, G] = fun(u);
N = numel(u);

disp('Testing the gradient');
for p = 2:12
    
    h = 10^(-p);
    Gh = zeros(N, 1);
    
    for i = 1:N
        uh = u; uh(i) = uh(i) + h;
        Fh = fun(uh);      
        Gh(i) = (Fh - F) / h;
    end
    
    disp([ 'p = ', num2str(p), ', err = ', num2str(norm(G- Gh, inf)) ]);
end

errG = G - Gh;

end
