

%% main choices of benchmark
% crack size
Kc = 2;
% data for toy-eam model
a = 4; b = 3; c = 20; rCutH = 1; 
% macro-deformation from ground state
F0 = [1, 0.02; 0, 1.02];
% directory name
dirname = 'bench_mcrack_2';


%% setup methods and method parameters
% reference atomistic region sizes (~powers of sqrt(2))
K0 = [3, 4, 6, 8, 11, 16, 22, 30];   
% K0 = [3, 4, 6, 8, 11];
% Full Atomistic with cut-off
meth(1).actype = 'atm';
meth(1).N = ceil(sqrt(2).^(6:12));
meth(1).fn = 'atm';
meth(1).alpha = 1.5;
% QCE Method
meth(2).actype = 'qce';
meth(2).K = K0 + 1;
meth(2).N = K0.^2 + 10;
meth(2).fn = 'qce';
meth(2).alpha = 1.5;
% B-QCE Method
meth(3).actype = 'bqc';
meth(3).blendtype = 'min-hessian';
meth(3).K = 2*K0;
meth(3).blendw = K0;
meth(3).N = K0.^2 + 10;
meth(3).fn = 'bqc_hess';
meth(3).alpha = 1.5;
% GRAC-23 Method
meth(4).actype = 'gqc23';
meth(4).K = K0+2;
meth(4).N = K0.^2 + 10;
meth(4).fn = 'gqc23';
meth(4).alpha = 1.5;
% REFERENCE
maxN = 0; for j = 1:4, maxN = max(maxN, max(meth(j).N)); end
% Option 1: GQC23
meth(5).actype = 'gqc23';
meth(5).K = 2*max(K0+2);
meth(5).N = max( (2*max(K0+2)).^2, 3 * maxN );
meth(5).fn = 'refsoln';
meth(5).alpha = 1.4;
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


%% create task list

ind = 0;
for jM = 1:length(meth)
  for jN = 1:length(meth(jM).N)
    ind = ind + 1;
    task(ind).actype = meth(jM).actype;
    task(ind).N = meth(jM).N(jN);
    task(ind).Kc = Kc;
    task(ind).model = model;
    task(ind).F0 = F0;
    task(ind).alpha = meth(jM).alpha;
    task(ind).jM = jM;
    task(ind).jN = jN;
    
    switch meth(jM).actype
      case 'atm'
        task(ind).K = task(ind).N;
      case 'qce'
        task(ind).actype = 'bqc';
        task(ind).blendtype = 'qce';
        task(ind).K = meth(jM).K(jN);
      case 'bqc'
        task(ind).blendtype = meth(jM).blendtype;
        task(ind).bw = meth(jM).blendw(jN);
        task(ind).K = meth(jM).K(jN);
      case 'gqc23'
        task(ind).K = meth(jM).K(jN);
      otherwise
        error('ERROR: unknown actype in bench_mcrack');
    end
    
    task(ind).filename = ['./out/', dirname, '/',  ...
                          meth(jM).fn, '_', ...
                          num2str(task(ind).N) ];
  end
end

ntasks = ind;

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
save(['out/', dirname, '/setup'], 'data', 'task', 'meth');

    
%% main benchmark loop
% create a storage function
parforsave = @(fn, geom, Y, E, task)(save(fn, 'geom', 'Y', 'E', 'task'));
% start loop
parfor t = 1:ntasks
  % display a message to show which task this is:
  disp([task(t).actype, ': N = ', num2str(task(t).N)]);
  % check whether file exists. 
  if exist([task(t).filename, '.mat'], 'file') ~= 2
    
    % create geometry
    geom = geom_2dtri_mcrack(task(t).Kc, task(t).K, task(t).N, rCutH, ...
      task(t).alpha);
    % prepare geom for AC-method
    switch task(t).actype
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
    disp(['   #DoFs : ', num2str(geom.nX)]);

    % setup optimisation
    addpath ./packages/popt
    opts = solve_opts();
    opts.actype = task(t).actype;
    opts.popt.disp = 0;
    % call the solver
    [Y, E] = solve_main(geom, task(t).model, task(t).F0, opts);
    opts.prec = 'none';
    opts.U0 = Y;
    [Y, E] = solve_main(geom, task(t).model, task(t).F0, opts);
    % save result
    parforsave(task(t).filename, geom, Y, E, task(t));
  end
end




