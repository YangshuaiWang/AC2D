% ex_mcrack
%
% most basic example to test codes on
%


%% Parameters
% type of ac-method:
actype = 'atm';       % {qce, bqc, atm, gqc23}   
% domain parameters
N = 50;         % approx. radius of domain (in atomic spacings)
Kc = 7;          % size of crack is 2*K0+1 atoms
% macroscopic strain relative to ground-state
F0 = [1, 0.02; 0, 1.03];
% model parameters: toy_eam_h
a = 4.0; b = 3; c = 10; rCutH = 1;

%% setup

% 1. Model
model = model_toyeam_h(a, b, c, rCutH);
model.bc = 'dir';

% now that we have defined the model, we can get the 
% macroscopic ground state > macro-scopic strain relative to Id
Fcb = get_cb_groundstate(model);
F0 = F0 * Fcb;

% 2. get some parameters that depend on the choice of the a/c method
switch actype
  case 'bqc'
    K0 = 4; bw = ceil(K0/2+1);
    K = K0 + bw + rCutH;
    blendtype = 'min-hessian';
  case 'qce'
    bw = 0;
    blendtype = 'qce';
    K = 6 + model.rCutH;
    actype = 'bqc';
  case 'gqc23'
    K = 6 + model.rCutH;
    if rCutH > 1
      error('ERROR: gqc23 requires rCutH = 1');
    end
  case 'atm'
    K = N;

  otherwise
    error('ERROR: unknown ac method');
end

% 3. generate geometry
geom = geom_2dtri_mcrack(Kc, K, N, rCutH, 1.3);
% actype-dependent geometry preparation
switch actype
  case 'bqc'
    geom = bqc_prep_geom(geom, 2, bw, blendtype);
    % geom = bqc_prep_geom(geom, 2, bw, blendtype);
  case 'gqc23'
    geom = gqc23_prep_geom(geom);
  case 'atm'
    geom.volX = ones(1, geom.nX);
end


%% solver
sopts = solve_opts();
sopts.visualize = false;
sopts.actype = actype;
sopts.popt.disp = 1;

profile on;
[Y, E] = solve_main(geom, model, F0, sopts);
profile off; profile report;


%% plot result
figure(1);
U = Y - F0 * geom.X;
nrmU = sqrt(sum(U.^2, 1));
% nrmDu = W1p_norm(geom.X, geom.T, U, inf);
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', U(1,:)');

figure(2);
if strcmp(actype, 'atm')
  plot(Y(1,:), Y(2,:), 'k.', 'Markersize', 20);
else
  trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', geom.volX);
end
view(2);



figure(3)
subplot(1,2,1);
nrmX = sqrt(sum(geom.X.^2,1));
r = sort(nrmX);
loglog(nrmX, nrmU, '.', r, max(nrmU) * r.^(-1), 'r-');
xlabel('|X|'); ylabel('|U|'); legend('|U|', 'r^{-1}');

subplot(1,2,2);
[rDu, Du] = get_grad_r(U, geom);
loglog(rDu, Du, '.', rDu, 0.02*max(Du) * rDu.^(-2), 'r-', ...
  rDu, 11*max(Du) * rDu.^(-2), 'r-');
xlabel('|X|'); ylabel('|DU|'); legend('|DU|', 'r^{-1}');

figure(4)
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', log(nrmU+exp(-10))');

%%

% figure(5);
% pub_vis_bqce(geom, Y, 7, 0, 0, 1.0);
% axis([-100, 100, -100, 100]);
% set(gca, 'XTick', [-100:50:100]);
% set(gca, 'YTick', [-100:50:100]);
% set(gca, 'Fontsize', 14);
% % xlabel('atomic units');
% % ylabel('atomic units');
% 
% figure(6);
% pub_vis_bqce(geom, Y, 100, 0, 0, 1);
% axis([0, 25, -12.5, 12.5]);
% set(gca, 'XTick', [0:5:25]);
% set(gca, 'YTick', [-10:5:10]);
% set(gca, 'Fontsize', 14);
% % xlabel('atomic units');
% % ylabel('atomic units');



