% A script to run the adaptive process based on A Posteriori Error
% Estimator and GQC23 Coupling.

% initialization
% ex_mcrack;

% parameter set
n = 0;
adapk = 1;
Tol = 1200; % 1500;
sample = 'GD'; % {'GD', 'ED'}
% dirname = sprintf('mcrackV%d-EAM', nV);
% dirname = sprintf('sample');
% dirname = sprintf('noDivFree');
dirname = 'bench_mcrack_initial';
load(['./out/', dirname, '/setup']);
load(meth(1).filename);
model = data.model;
model.Kc = data.Kc;
iRx = find(geom.X(2,:) == 0 & geom.X(1,:) > 0);
Rx = geom.X(1,iRx);
geom.iRx = iRx;
geom.Rx = Rx;
geom.adapk = adapk;
adapflow(1).refX = [];
adapflow(1).Ixh = 0;
fn = sprintf('out/%s/%s', dirname, sample);
while n ~= 1
    [geom, model, Y, adapflow] = solve_adaptive_main_ys(geom, model, Y, dirname, sample, fn, adapflow);
    geom.iRx = iRx;
    geom.Rx = Rx;
%     keyboard;
    dof = geom.nX;
    adapk = adapk + 1;
    geom.adapk = adapk;
    if dof >= Tol
        n = n + 1;
    end
end
