% function error = error_main(geom, square, STh, STe)
% Program: compute error estimate.
% Input:
%     geom: geometry structure.
%     square: \epsilon scale geometry structure.
%     STh: stress tensor of T_h
%     STe: stress tensor of T_epsilon
% Output: 
%     err: totoal error. 
% Author : M.Liao
% Version:
%     1.1 correct those constant expect stability constant; apply area to 
%       compute coarsen error instead of max h. 11-Nov-2014.
%     1.11, type added to decide what error to be sort 21-Mar-2015.
%     1.2 applying lambda and compute prop     Mar-16-2016
function error = error_main(geom, square, STh, STe)
% This file is used to compute errors 

% check input
% We might better leave these work to the partitional files.
% TODO: 1. test routine, 2. we might check U, as well. 
if nargin == 0
  test_error_main();
  return;
end
% We should emphaizt that the values of epsilon and h which initialized
% as following.

% parameter set
err_mdl = error_model(geom, square, STh, STe);
sum_mdl = sum(err_mdl);

if sum_mdl ~= 0
    err_mdl = sqrt(sqrt(3)/4)/sqrt(sum_mdl)*err_mdl;
end
% compute coarsening error
err_crsn = error_coarsen(geom, STh);
sum_crsn = sum(err_crsn);

if sum_crsn ~= 0
    err_crsn = sqrt(3)/sqrt(sum_crsn)*err_crsn;
end
% ====================================================================== %
% lambda applied
lambda = sum(err_crsn)/sum(err_mdl);
% if lambda >= 1;
%     err_mdl = lambda*err_mdl;
% else
%     err_crsn = err_crsn/lambda;
% end
% lambda = 16.6/sqrt(3);
% err_mdl = lambda*err_mdl;
% err_mdl = 18.563*err_mdl;
% ====================================================================== %
% compute & sort error
err = err_mdl + err_crsn;

% ====================================================================== %
% apply a constant to feature the elements to refine!
% [err_rfn idx] = sort(err, 'descend');
% err_total = cumsum(err_rfn);
% T_rfn = find(err_total > err_total(end)*0.5, 1);
% T_rfn = idx(1:T_rfn);
% ====================================================================== %
% apply mean value of nonzero estimators to feature to-be refined elements!
errTmp = err;
errTmp(err == 0) = [];
Tol = mean(errTmp);
T_rfn = find(err > Tol);
% ====================================================================== %

iX = find(geom.di == 0);
iT_intf = geom_2dtri_aIXT(geom, iX);
iT = intersect(iT_intf, T_rfn);
prop = sum(err(iT))/sum(err);

error.err = err;
error.total = sum(err);
error.model = err_mdl;
error.coarsen = err_crsn;
% error.T_order = idx;
error.T_refine = T_rfn;
error.lambda = lambda;
error.prop = prop;

end


%% TEST ROUTINE
function test_error_main()
geom = geom_2dtri_mcrack(3, 3, 20, 1, 1.5);
geom = geom_analyze(geom, 1);
geom = gqc23_prep_geom(geom);
geom = geom_2dtri_h(geom);
model = model_toyeam_h(4, 3, 0, 1);%6*exp(-3), 1);

% geom = geom_2dtri_mcrack(1, 3, 10, 1, 1.5);
% % geom = bqc_prep_geom(geom, 1, 2, 'min-hessian');
% geom = gqc23_prep_geom(geom);
% model = model_toyeam_h(4, 3, 0, 6*exp(-3), 1);
geom.plot(geom);
% geom = geom_2dtri_h(geom);
for i = 1:geom.nT
    text(geom.heart(1,i), geom.heart(2,i), num2str(i));
end
A = [1 .001; .002 1];
U = A*geom.X + 1e-2*rand(size(geom.X));

model.F0 = A;
model.rho0 = 6*exp(-3);
model.W0 = 0.0;

square = geom_2dtri_scross(geom);
square = geom_interp_scross(geom, square, U, model);
square = geom_scross_neigs(geom, square, U, model);
STh = get_stress_tensor(geom, model, U);
% [~, STe] = scross_stress_tensor(square, model, A);
[~, STe] = scross_stress_tensor(square, model, geom, U);
[S, ufs, divFC] = error_div_free(geom, square, STh, STe);
STh = gqc23_stress_tensor(STh, S, ufs, divFC);
error = error_main(geom, square, STh, STe);
disp(num2str(error.total));
end
