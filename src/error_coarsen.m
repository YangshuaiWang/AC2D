% function err_crsn = get_coarsen_error(geom, sigma)
% 
% input: 
%     geom, geometry structure.
%     sigma, stress tensor.
% output: 
%     err_crsn, coarsening error.
% Version: 
%     1.1 Iteration over edges instead of triangles. remove function
%     get_neig_tris(geom);

function err_crsn = error_coarsen(geom, sigma)
% test routine
if nargin == 0
  test_get_coarsen_error();
  return;
end

% check input
if (~isfield(geom, 'volT'))
  error('ERROR: must call gqc23_prep_geom before get_model_error!');
end
%% main iteration
err_crsn = zeros(1, geom.nT);
E = geom.E;
volT = geom.volT;
aET = geom.aET;
for i = 1:length(E)
%     if i == 459; keyboard; end
    T = aET(:,i);
%     if any(T == [0,0]') || any(volT(T)==0);
%         continue;
%     end
%     if any(volT(T)<0.5+eps); continue; end;

    if any(T == [0,0]') || all(volT(T)==0);
        continue;
    end
    e = geom.X(:,E(:,i));
    e = diff(e,[],2);
    hf = norm(e);
% ====================================================================== %
% older version
%     dw = reshape(sigma(:,:,T),4,2);
%     dw = diff(dw, [], 2);
%     dw_2 = sum(dw.^2);
% 
%     if any(volT(T)<0.5)
%         volT(T(volT(T)<0.5)) = 0;
%     end
%     err_crsn(T) = err_crsn(T) + dw_2*volT(T);
% 
% ====================================================================== %
    dw = sigma(:,:,T);
    jump = diff(dw,[],3)*null(e');
    jump = sum(jump.^2);
    err_crsn(T) = err_crsn(T) + .5*hf*jump;
end
err_crsn(volT == 0) = 0;
end
%% Test routine
function test_get_coarsen_error()
geom = geom_2dtri_mcrack(3, 3, 8, 1, 1.5);
% geom = bqc_prep_geom(geom, 1, 2, 'min-hessian');
geom = gqc23_prep_geom(geom);
% geom = geom_2dtri_mapping(geom);
model = model_toyeam_h(4, 3, 0, 6*exp(-3), 1);
geom.plot(geom);
% for i = 1:geom.nT
%     text(geom.heart(1,i), geom.heart(2,i), num2str(i));
% end
A = [1 .001; .002 1];
U = A*geom.X;
tic;
Sigma = get_stress_tensor(geom, model, U);
toc; tic;
err_mdl = error_coarsen(geom, Sigma);
toc;
str = sprintf('Test finished: maxinum error = %3.5f', max(err_mdl));
disp(str);
end
