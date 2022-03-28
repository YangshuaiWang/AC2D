
clear; clc

model = model_toyeam_h(4.0, 2.0, 10, 1);
N = 30;
K = 8;
F0 = [1, 0; 0, 1];
geom = geom_2dtri_hexagon_edge(N, K, 2, 'dir');
% geom = geom_2dtri_hexagon(N, K, 2, 'dir');
% geom = geom_create_vacancies(geom, 1);
geom = geom_analyze(geom);
% geomb = geom;
geom = gqc23_prep_geom(geom);

model.bc = 'dir';
Fcb = get_cb_groundstate(model);
F0 = F0 * Fcb;
model.F0 = F0;

% setup optimisation
addpath ./packages/popt
opts = solve_opts();
opts.actype = 'gqc23';
opts.popt.disp = 0;
% call the solver
[Y, E] = solve_main(geom, model, model.F0, opts);
opts.prec = 'none';
opts.U0 = Y;
[Y, E] = solve_main(geom, model, model.F0, opts);
U = Y - model.F0*geom.X;
trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', abs(U(1,:)));


% geomb = bqc_prep_geom(geomb, 1, 1, 'linearH');
% 
% geom.plot(geom);
% disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);
% for n = 1:geom.nX
%   text(geom.X(1,n)+0.1, geom.X(2,n), num2str(geom.volX(n)));
% end


% % make the ghost-force-test:
% U = (eye(2) + 0.02 * rand(2)) * geom.X;
% model.bc = 'dir';
% [E, dE] = gqc23_energy(U, geom, model); 
% % energy-test
% Eb = bqc_energy(U, geomb, model);
% disp(['Energy-Test: |Eq - Eb| / |Eb| = ', num2str(abs(E - Eb) / abs(Eb))]);
% % ghost-force test
% dE = reshape(dE, model.rDim, geom.nX); 
% gf = dE;
% gf(:, 1:28) = 0;
% dE = dE - gf;
% dE(:, geom.iBdry) = 0;
% nrmG = norm(dE(:), inf);
% disp(['Ghost force test: ||dE(yA)||_inf = ', num2str(nrmG)]);
% trisurf(geom.T', geom.X(1,:)', geom.X(2,:)', abs(dE(1,:)));
% 
% % energy functional
% model.bc = '';
% fcnl = @(U_)(gqc23_energy(U_, geom, model));

% % base point
% U = geom.X + 0.01 * rand(model.rDim, geom.nX);

% % call finite difference test
% addpath ./packages/popt
% disp('--------------------------------------------');
% disp('   Testing dE assembly in gqc23_energy.m');
% disp('--------------------------------------------');
% test_derivatives(fcnl, U, 1);
% disp('--------------------------------------------');

%
% geom.volX = ones(1, geom.nX);
% geom.plot(geom)

% Fcb = get_cb_groundstate(model);
% F0 = F0 * Fcb;
% model.F0 = F0;
% 
% % setup optimisation
% addpath ./packages/popt
% opts = solve_opts();
% opts.actype = 'atm';
% opts.popt.disp = 0;
% % call the solver
% [Y, E] = solve_main(geom, model, model.F0, opts);
% opts.prec = 'none';
% opts.U0 = Y;
% [Y, E] = solve_main(geom, model, model.F0, opts);

% function [E, dE] = test_edge_solve(U, geom, model)
% 
% if nargin == 0
%   test_atm_energy_edge();
%   return;
% end
% 
% end

% atm prep_geom
% geom.volX = ones(1, geom.nX);
% F0 = [1, 0; 0, 1];



% [Y, E] = solve_main(geom, model, F0, opts)

% function test_atm_energy_edge()
% % define model and geometry
% 
% model = model_toyeam_h(4.0, 2.0, 10, 1);
% N = 30;
% K = 8;
% geom = geom_2dtri_hexagon_edge(N, K, 2, 'dir');
% geom = geom_analyze(geom);
% geom.volX = ones(1, geom.nX);
% model.F0 = [1, 0; 0, 1];
% 
% disp(['nX = ', num2str(geom.nX), '; nT = ', num2str(geom.nT)]);
% 
% % wrapper functional
%   function [E, dE] = atm_wrapper(W, geom, model, iFree)
%     Y = geom.X;
%     Y(:, iFree) = reshape(W, 2, length(iFree));
%     if nargout == 1
%       E = atm_energy(Y, geom, model);
%     else
%       [E, dE] = atm_energy(Y, geom, model);
%       dE = reshape(dE, 2, geom.nX);
%       dE = dE(:, iFree);
%       dE = dE(:);
%     end
%   end
% 
% % energy functional
% % iFree = setdiff(1:geom.nX, geom.iPer(1,:));
% if isfield(geom, 'iDir')
%   iFree = setdiff(1:geom.nX, geom.iDir);
% else
%   iFree = setdiff(1:geom.nX, geom.iBdry);
% end
% 
% fcnl = @(W_)(atm_wrapper(W_, geom, model, iFree));
% 
% % base point
% U = geom.X + 0.01 * rand(model.rDim, geom.nX);
% W0 = U(:, iFree);
% 
% % finite difference test
% addpath ./packages/popt
% disp('--------------------------------------------');
% disp('   Testing dE assembly in assemble_a.m');
% disp('--------------------------------------------');
% test_derivatives(fcnl, W0(:), 1);
% disp('--------------------------------------------');
% 
% end
