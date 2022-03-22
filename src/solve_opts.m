% opts = solve_opts
%
% TODO: write this routine properly
%
% opts.actype  {qce, bqc, gqc23}
% opts.prec    {laplace, none}
% opts.solver  {popt}
%     opts.popt : popt options
%

function opts = solve_opts()

opts.actype = 'bqc';
opts.prec = 'laplace';
opts.solver = 'popt';
opts.visualize = false;
opts.U0 = [];
opts.post_fsolve = false;

addpath ./packages/popt/
opts.popt = poptls_opts('pr+', 'rough');
opts.popt.g_tol = 1e-5; 
opts.popt.g_tol_inf = 1e-6;
opts.popt.x_tol = 1e-5; 
opts.popt.x_tol_inf = 1e-7;
opts.popt.disp = 0;
end
