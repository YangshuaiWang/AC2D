% function err_mdl = get_model_error(geom, square, STh, STe)
% Input: geom, geometry structure;
%        modle, energy model structure;
%        STh, sigma^{h} by calling get_stress_tensor;
%		 square, structure save info for T_e;
%		 STe, sigma^{a} by calling scross_stress_tensor.
% Output: err_mdl: modelling error.
% Author: M.Liao
% Version: 1.1 square invoked to make the process clear. Nov-28-2014.
%          1.11 assemble to T by half. Jan-26-2015.
%          1.12 assemble to T by wight as |T|/w(T_e). Actually not.. by 1..
%          1.2 assemble according to modified estimator. i.e. eta_1.
%          1.21 replace sigma^{a} by STe. 	Oct-19-2015.
function err_Th = error_model(geom, square, STh, STe)
% compute modelling error. By geometry structure and displacement.

% TODO: test routine.
if nargin == 0
  test_get_model_error(); % Does not work now.. 
  return;
end
% check input

% wheather geom_analyze been called checked in bqc_prep_geom.
if (~isfield(geom, 'volT'))
  error('ERROR: must call bqc_prep_geom before get_model_error!');
end

err_Th = zeros(1,geom.nT);
T = square.T;
OA = square.oa;
aIoa = square.aIoa;  
aTI = square.aTI;
TI = square.iTI;
TeI = square.iTeI;
Tec = setdiff(1:size(T,2), TeI);
% main loop over small triangles to compute difference of stress and
% assemble to overlapping elements. p.s. small triangles are divided into 
% Tec and Teic as before, for Teic represents /mathcal(T)^{i,c}, with Tec 
% its complement.
for i = 1:numel(TeI)
	se = STe(:,:,TeI(i));
	sh = STh(:,:,TI(i));
	ds = reshape((se-sh),4,1);
	ds = sum(ds.^2);	
	err_Th(TI(i)) = err_Th(TI(i)) + ds;
end

for i = 1:numel(Tec)
%     t = T(:,i);
    i_sqr = aTI(1,i);
    psn = aTI(2,i);
    oa = OA(:,aIoa == i_sqr);
    Th_idx = oa(1,:);
    wght = oa(psn,:)./sum(oa(psn,:));
    if numel(wght) == 1; continue; end;
    if any(abs(wght)<eps)||any(isnan(wght)); continue; end;
    
    sh = STh(:,:,Th_idx);
    lnth = size(sh,3);
    
    sh = reshape(sh,4,lnth);
	se = STe(:,:,i);
    se = reshape(se,4,1);
% ====================================================================== %    
%     computation of eta, for checking.
    ds = se - sum(kron(wght, ones(4,1)).*sh, 2);
    ds = ds.^2;
    ds = sum(ds);
%     if find(Th_idx == 35, 1); keyboard; end
% assembling
if lnth > 1;
%     pair = combntns(1:lnth,2);
%     dS = kron(wght, ones(4,1)).*(repmat(se,1,lnth) - sh);
%     sum_dS = sum(dS.^2);
%     err_Th(Th_idx) = err_Th(Th_idx) + sum_dS;
%     for j = 1:size(pair,1);
%         p = pair(j,:);
%         a = dS(:,p(1));
%         b = dS(:,p(2));
%         err_Th(Th_idx(p(1))) = err_Th(Th_idx(p(1))) + ...
%             sum_dS(p(1))/sum(sum_dS(p))*sum(2*(a.*b));
%         err_Th(Th_idx(p(2))) = err_Th(Th_idx(p(2))) + ...
%             sum_dS(p(2))/sum(sum_dS(p))*sum(2*(a.*b));
%     end
    dS = kron(wght, ones(4,1)).*(repmat(se,1,lnth) - sh);
    sum_dS = sum(dS.^2);
    err_Th(Th_idx) = sum_dS/sum(sum_dS)*ds;
    continue
end
% ====================================================================== %
    err_Th(Th_idx) = err_Th(Th_idx) + ds;
end
end

%% test routine
function test_get_model_error()
geom = geom_2dtri_mcrack(3, 3, 20, 1, 1.5);
% geom = bqc_prep_geom(geom, 1, 2, 'min-hessian');
model = model_toyeam_h(4, 3, 0, 6*exp(-3), 1);
% geom = geom_atomistic_expansion(geom, 2);
geom = geom_analyze(geom, 1);
geom = gqc23_prep_geom(geom);
geom.plot(geom);
for i = 1:geom.nT
    t = geom.T(:,i); x = sum(geom.X(1,t))/3; y = sum(geom.X(2,t))/3;
    text(x, y, num2str(i));
end
A = [1 .001; .002 1];
U = A*geom.X;
tic;

STh = get_stress_tensor(geom, model, U);
square = geom_2dtri_scross(geom);
square = geom_interp_scross(geom, square, U, model);
toc; tic;
err_Th = error_model(geom, model, STh, square);
toc;
str = sprintf('Test finished: maxinum error = %3.5f', max(err_Th));
disp(str);
end
