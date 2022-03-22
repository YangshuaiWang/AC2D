

% dirname = 'bench_mcrack_4'    %#ok
dirname = 'bench_mcrack_4_ys'

% load setup data
load(['out/', dirname, '/setup']);

%% Compute errors


% If the variable "errors_computed" does not exist, then we first 
% need to compute them
if ~exist('errors_computed', 'var')
  
  % load reference solution and store required reference values
  load(task(end).filename, 'geom', 'Y', 'E');
  U_ref = Y - data.F0 * geom.X;
  geom_ref = geom;
  E_ref = E;
  ref_nrm2 = W1p_norm(geom.X, geom.T, U_ref, 2);
  ref_nrminf = W1p_norm(geom.X, geom.T, U_ref, inf);
  ref_nrm1 = W1p_norm(geom.X, geom.T, U_ref, 1);
  
  % loop through computations
  for t = 1:(length(task)-1)
    disp(task(t).filename);
    % load solution
    load(task(t).filename, 'geom', 'Y', 'E');
    U = Y - data.F0 * geom.X;
    % E = E - geom_volume(geom) * W0;
    
    % extract number of dofs
    meth(task(t).jM).dof(task(t).jN) = geom.nX;
    
    % W1p-errors
    [err1, err2, errinf] = err_W1p(U, geom, U_ref, geom_ref);
    meth(task(t).jM).Err1(task(t).jN) = err1 / ref_nrm1;
    meth(task(t).jM).Err2(task(t).jN) = err2 / ref_nrm2;
    meth(task(t).jM).ErrInf(task(t).jN) = errinf / ref_nrminf;

    % energy-error: !!!TODO!!!
    data.meth(task(t).jM).ErrE(task(t).jN) = 1;
  end

  errors_computed = true;
  save(['out/', dirname, '/setup'], 'data', 'task', ...
    'meth', 'errors_computed');
end


% TODO:
% % % compute the energy of the ground-state to compare against
% % volDom = 6 * sqrt(3)/2 * data.N^2;
% % W0 = data.model.Wfun(model, data.F0);
% % E_ref = E_ref - volDom * W0;
% 
% % subtract the energy of the reference state
% W0 = data.model.Wfun(data.model, data.model.F0);
% E_ref = E_ref - geom_volume(geom_ref) * W0;





%% Plots  
mrkrs = {'o-', 's-'};
cols = {[0,0,0], [1,0,0]};
meth(1).actype = 'ATM';
meth(2).actype = 'GRAC';
dof = [1e2, 1e5];
lw = 2.0; lwth = 1.0; ms = 12;
mrk_rate = 'k:';


nac = length(meth)-1;

% W12-error
figure(1); 
for j = 1:nac
  loglog(meth(j).dof, meth(j).Err2, mrkrs{j}, ...
    'Linewidth', lw, 'Markersize', ms, 'Color', cols{j}); hold on;
end
loglog(dof, 0.35*(dof/dof(1)).^(-1/2), 'k:', 'Linewidth', lwth);

currfig = gcf;

annotation(currfig,'textbox',...
  [0.637363636363636 0.68 0.3 0.0389016018306637],...
  'String',{'(# DoFs)^{-1/2}'},...
  'FitBoxToText','off', 'Linestyle', 'none', 'Fontsize', 20);

loglog(dof, 0.3*(dof/dof(1)).^(-1), 'k:', 'Linewidth', lwth); 
annotation(currfig,'textbox',...
  [0.68 0.34 0.3 0.0389016018306637],...
  'String',{'(# DoFs)^{-1}'},...
  'FitBoxToText','off', 'Linestyle', 'none', 'Fontsize', 20);

loglog([1, 1e6], [1e-2, 1e-2], 'm-'); 


hold off;
xlabel('# DoFs', 'Fontsize', 18); 
ylabel('$| u_{\mathrm{ac}} - u |_{1,2} \, / \, | u |_{1,2}$', ...
  'Interpreter', 'latex', 'Fontsize', 18);
% axis([dof(1) / 1.5, dof(2) * 1.5, 5e-3, 5e-1]);
axis([dof(1) / 1.5, dof(2) * 1.5, 5e-3, 1]);
% legend('QCE', 'Linear Blending', 'Smooth Blending', 
% 'DoF^{-1/2}', 'Location', 'Southwest');
set(gca, 'XMinortick', 'off', 'YMinorTick', 'off', 'Fontsize', 18);

for j = 1:nac
  acmethods{j} = meth(j).actype; %#ok
end
legend(acmethods, 'Location', 'Southwest');

% W1inf-error
figure(2);
for j = 1:nac
  loglog(meth(j).dof, meth(j).ErrInf, mrkrs{j}, ...
    'Linewidth', lw, 'Markersize', ms, 'Color', cols{j}); hold on;
end
% loglog(dof, (dof/dof(1)).^(-1), 'k:', 'Linewidth', lwth); 
% loglog(dof, (dof/dof(1)).^(-3/2), 'k:', 'Linewidth', lwth); 
hold off;
xlabel('# DoFs', 'Fontsize', 18); 
ylabel('$| u_{\mathrm{ac}} - u |_{1,\infty} \, / \, | u |_{1,\infty}$', ...
  'Interpreter', 'latex', 'Fontsize', 18);
% axis([dof(1) / 1.5, dof(2) * 1.5, 2e-3, 2e-1]);
axis([dof(1) / 1.5, 1e4, 4e-3, 1.5e-1]);
% legend('QCE', 'Linear Blending', 'Smooth Blending', 'DoF^{-1/2}', 'Location', 'Southwest');
set(gca, 'XMinortick', 'off', 'YMinorTick', 'off', 'Fontsize', 18);
legend(acmethods, 'Location', 'Southwest');


% % W11-error
% subplot(2,2,3);
% for j = 1:nac
%   loglog(meth(j).dof, meth(j).Err1, mrkrs{j}, ...
%     'Linewidth', lw, 'Markersize', ms); hold on;
% end
% loglog(dof, (dof/dof(1)).^(-1/4), 'k:', 'Linewidth', lwth); 
% loglog(dof, (dof/dof(1)).^(-1/2), 'k:', 'Linewidth', lwth); 
% hold off;
% xlabel('# DoFs'); 
% ylabel('$| u_{\mathrm{ac}} - u |_{1,1} \, / \, | u |_{1,1}$', ...
%   'Interpreter', 'latex');
% axis([dof(1) / 1.5, dof(2) * 1.5, 1e-1, 1]);
% % legend('QCE', 'Linear Blending', 'Smooth Blending', 'DoF^{-1/2}', 'Location', 'Southwest');
% set(gca, 'XMinortick', 'off', 'YMinorTick', 'off', 'Fontsize', 14);

return;

%% W12-error, 4 separate figures
close all;
nac_ = nac;
acmethods = [];
for nac = 1:nac_

  
% W12-error
currfig = figure(nac);
for j = 1:nac
  loglog(meth(j).dof, meth(j).Err2, mrkrs{j}, ...
    'Linewidth', lw, 'Markersize', ms, 'Color', cols{j} );  hold on;
end

loglog([1, 1e6], [1e-2, 1e-2], 'm-', 'Linewidth', 0.5); 
% loglog([1, 1e6], [1e-1, 1e-1], 'm-', 'Linewidth', 0.5);
% loglog([1, 1e6], [3e-2, 3e-2], 'm-', 'Linewidth', 0.5); 

loglog(dof, 0.35*(dof/dof(1)).^(-1/2), 'k:', 'Linewidth', lwth);
annotation(currfig,'textbox',...
  [0.637363636363636 0.68 0.3 0.0389016018306637],...
  'String',{'(# DoFs)^{-1/2}'},...
  'FitBoxToText','off', 'Linestyle', 'none', 'Fontsize', 20);

if nac >= 4
loglog(dof, 0.3*(dof/dof(1)).^(-1), 'k:', 'Linewidth', lwth); 
annotation(currfig,'textbox',...
  [0.68 0.34 0.3 0.0389016018306637],...
  'String',{'(# DoFs)^{-1}'},...
  'FitBoxToText','off', 'Linestyle', 'none', 'Fontsize', 20);
end


hold off;
xlabel('# DoFs', 'Fontsize', 20); 
ylabel('$| u_{\mathrm{ac}} - u |_{1,2} \, / \, | u |_{1,2}$', ...
  'Interpreter', 'latex', 'Fontsize', 20);
axis([dof(1) / 1.5, dof(2) * 1.5, 4e-3, 4e-1]);
% legend('QCE', 'Linear Blending', 'Smooth Blending', 
% 'DoF^{-1/2}', 'Location', 'Southwest');
set(gca, 'XMinortick', 'off', 'YMinorTick', 'off', 'Fontsize', 18);

for j = 1:nac
  acmethods{j} = meth(j).actype; %#ok
end
legend(acmethods, 'Location', 'SouthWest');

end


% % energy-error
% subplot(2,2,4);
% for j = 1:nac
%   loglog(meth(j).dof, meth(j).ErrE, mrkrs{j}, ...
%     'Linewidth', lw, 'Markersize', ms); hold on;
% end
% loglog(dof, (dof/dof(1)).^(-1), 'k:', 'Linewidth', lwth); 
% loglog(dof, (dof/dof(1)).^(-2), 'k:', 'Linewidth', lwth); 
% hold off;
% xlabel('# DoFs'); 
% ylabel('$| u_{\mathrm{ac}} - u |_{1,1} \, / \, | u |_{1,2}$', ...
%   'Interpreter', 'latex');
% axis([dof(1) / 1.5, dof(2) * 1.5, 1e-2, 1]);
% % legend('QCE', 'Linear Blending', 'Smooth Blending', 'DoF^{-1/2}', 'Location', 'Southwest');
% set(gca, 'XMinortick', 'off', 'YMinorTick', 'off', 'Fontsize', 14);
% 


