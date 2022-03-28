
% 2022/03 test

%% main choices of benchmark
% crack size
Kc = 3;

% whole domain size
N = 100;

% data for toy-eam model
a = 4; b = 3; c = 20; rCutH = 1; 
% macro-deformation from ground state
F0 = [1, 0.02; 0, 1.02];
% directory name
dirname = 'bench_mcrack_initial';

%% setup methods and method parameters
% reference atomistic region sizes (~powers of sqrt(2))
% K0 = [3, 4, 6, 8, 11, 16, 22, 30];   
K0 = 5; % [3, 4, 6, 8, 11];
% GQC23
meth(1).actype = 'gqc23';
meth(1).K = K0;
meth(1).N = N;
meth(1).fn = 'gqc23';
meth(1).alpha = 1.5;
meth(1).filename = ['./out/', dirname, '/', meth(1).fn];
% REFERENCE
maxN = 0; maxN = max(maxN, max(meth(1).N));
% Option 1: GQC23
meth(2).actype = 'gqc23';
meth(2).K = 2*max(K0+15);
meth(2).N = max( (2*max(K0+2)).^2, 3 * maxN );
meth(2).fn = 'refsoln';
meth(2).alpha = 1.4;
meth(2).filename = ['./out/', dirname, '/', meth(2).fn];
% % Option 2: ATM
% meth(5).actype = 'atm';
% meth(5).N = ceil(1.5 * maxN);
% meth(5).fn = 'refsoln';


%% more setup
% model
model = model_toyeam_h(a, b, c, rCutH);
model.bc = 'dir';
max_dh = model.rCutH;
% CB reference state to modify F0
Fcb = get_cb_groundstate(model);
F0 = F0 * Fcb;
model.F0 = F0;

%% save setup
% global data
data.bc = 'dir';
data.F0 = F0;
data.a = a; data.b = b; data.c = c; data.rCutH = rCutH;
data.K0 = K0;
data.Kc = Kc;
data.model = model;
data.dirname = dirname;

% save both data and task in the setup.mat
mkdir(['out/', dirname]);
save(['out/', dirname, '/setup'], 'data', 'meth');


%% main benchmark loop
% create a storage function
parforsave = @(fn, geom, Y, E, meth)(save(fn, 'geom', 'Y', 'E', 'meth'));
% start loop
for t = 1:2
  % display a message to show which task this is:
%   disp([task(t).actype, ': N = ', num2str(task(t).N)]);
    fprintf('\t%s\n', meth(t).fn);
  % check whether file exists. 
%   if exist([task(t).filename, '.mat'], 'file') ~= 2
    
    % create geometry
    geom = geom_2dtri_mcrack(data.Kc, meth(t).K, meth(t).N, rCutH, ...
      meth(t).alpha);
    % prepare geom for AC-method
    switch meth(t).actype
      case 'atm'
        geom.volX = ones(1, geom.nX);
      case 'qce'
        geom = bqc_prep_geom(geom, rCutH, 1, 'qce');
      case 'bqc'
        geom = bqc_prep_geom(geom, rCutH, task(t).bw, task(t).blendtype);
      case 'gqc23'
        geom = gqc23_prep_geom(geom);
      otherwise
        error('ERROR: unknown actype in bench_vac');
    end
    
    % output number of DOFs
%     disp(['   #DoFs : ', num2str(geom.nX)]);
    fprintf('\t#DoFs  : %d\n', geom.nX);

    % setup optimisation
    addpath ./packages/popt
    opts = solve_opts();
    opts.actype = meth(t).actype;
    opts.popt.disp = 0;
    % call the solver
    [Y, E] = solve_main(geom, data.model, data.F0, opts);
    opts.prec = 'none';
    opts.U0 = Y;
    [Y, E] = solve_main(geom, data.model, data.F0, opts);
    % save result
    parforsave(meth(t).filename, geom, Y, E, meth(t));
 end
% end

% %% create task list
% 
% ind = 0;
% for jM = 1:length(meth)
%   for jN = 1:length(meth(jM).N)
%     ind = ind + 1;
%     task(ind).actype = meth(jM).actype;
%     task(ind).N = meth(jM).N(jN);
%     task(ind).Kc = Kc;
%     task(ind).model = model;
%     task(ind).F0 = F0;
%     task(ind).alpha = meth(jM).alpha;
%     task(ind).jM = jM;
%     task(ind).jN = jN;
%     
%     switch meth(jM).actype
%       case 'atm'
%         task(ind).K = task(ind).N;
%       case 'qce'
%         task(ind).actype = 'bqc';
%         task(ind).blendtype = 'qce';
%         task(ind).K = meth(jM).K(jN);
%       case 'bqc'
%         task(ind).blendtype = meth(jM).blendtype;
%         task(ind).bw = meth(jM).blendw(jN);
%         task(ind).K = meth(jM).K(jN);
%       case 'gqc23'
%         task(ind).K = meth(jM).K(jN);
%       otherwise
%         error('ERROR: unknown actype in bench_mcrack');
%     end
%     
%     task(ind).filename = ['./out/', dirname, '/',  ...
%                           meth(jM).fn, '_', ...
%                           num2str(task(ind).N) ];
%   end
% end
% 
% ntasks = ind;



    





