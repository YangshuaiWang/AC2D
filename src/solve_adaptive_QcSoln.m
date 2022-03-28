% function [Y, Err, geom] = ex_mcrack_posteriori(geom, K2)
% function to apply a posteriori estimator
% Input:
%     geom, proper geometry structure.
%     K2, number of atom layers to be added.
% Output:
%     Y, displacements solved
%     Err, error computed invoking error_main
%     geom, new proper geometry structure.
% function [Y, Err, geom, T_rfn] = ex_mcrack_posteriori(geom, model, K2)
function [Y, E, geom] = solve_adaptive_QcSoln(geom, model, F0, actype)
% F0 = [1.0, 0.03; 0, 1.03];
% Fcb = get_cb_groundstate(model);
% F0 = F0 * Fcb;
% actype = 'gqc23';

% geom = geom_atomistic_expansion(geom, K2);
% geom = geom_2dtri_mapping(geom);
if nargin == 3;
    actype = 'gqc23';
end

geom = geom_analyze(geom, 1);

if strcmp(actype, 'gqc23')
    if isfield(geom, 'volX')
        geom = rmfield(geom, {'volX','volT'});
    end
    geom = gqc23_prep_geom(geom);
end
if strcmp(actype, 'atm')
    geom.volX = ones(1, geom.nX);
    geom.volX(geom.di <= model.rCutH-1) = 0.0;
end

sopts = solve_opts();
sopts.visualize = false;
sopts.actype = actype;
sopts.popt.disp = 0;
% sopts.post_newton = true;
% sopts.popt.Psolve = 'cg';
% sopts.popt.nCG_P = 30;

disp('Nonlinear solver...'); tic;
[Y, E] = solve_main(geom, model, F0, sopts); toc

% [Err T_rfn] = error_main(geom, model, Y);
  
end