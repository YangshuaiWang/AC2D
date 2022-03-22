% [x, info] = popttr0(ffun, Pfun, x0, options)
% Preconditioned trust region method preliminary implementation
%
% Input:
%   ffun     : objective function; [f, g] = ffun(x)
%   Pfun     : preconditioner; P = Pfun(x), or constant matrix P, or []
%   x0       : initial guess
%   options  : structure (use poptls_options for now)
%
% Output:
%   x    : final iterate, empty array if psd was unsuccesful
%   info : structure with the following content
%            .n_it : number of iterations
%            .n_f  : number of function evaluations
%            .n_p  : number of preconditioner applications
%            .errmsg
%            .n_reject


%          if opts.hist = true, the following history is stored in info:
%            .h_res, .h_res_inf, .h_step, .h_step_inf
%            .h_f, .h_n_f, .h_n_it


function [x, info] = popttr0(ffun, Pfun, x, opts)


%% *********************  Initialization  *************************

% initialize info structure
info.errno = 0; 
info.errmsg = 'success';
info.n_it = 0; 
info.n_f = 0; 
info.n_p = 0;
info.n_reject = 0;
% first function evaluation
[f, g] = ffun(x);
info.n_f = info.n_f + 1;
% precondition the gradient
if opts.apply_P
  gP = Pfun(x, g);
else
  % compute preconditioner
  P = eval_P(Pfun, x, opts);
  % precondition the first gradient
  gP = apply_P(P, g, opts);
end
info.n_p = info.n_p + 1;
% auxiliary scalar: g'*gP = g'*P^{-1}*g = gP'*P*gP
g_dot_gP = g'*gP;
% initialize trust region radius
DELTA = min(1, sqrt(abs(gP'*g)));

%% ********************* Main Loop ******************
res = sqrt(g_dot_gP);
res_inf = norm(g, inf);
step = 0;
step_inf = 0;
f_old = f + opts.f_tol + 1;
if opts.disp >= 1
    fprintf(['-----------------------------------------------------' ...
             '--------------------------------\n']);
    fprintf('%5s %6s %11s %9s %9s %9s %9s %9s %9s\n', ...
      'n_it', 'n_f', 'f', 'Res(P)', 'Res(inf)', 'Step(P)', 'Step(Inf)', ...
      'DELTA', 'rho');
    fprintf(['-----------------------------------------------------' ...
             '--------------------------------\n']);
    fprintf('%5d %6d %11.3e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n', ...
      info.n_it, info.n_f, f, res, res_inf, step, step_inf, DELTA, 0);
end

% Notes on the termination criterion:
%   res = sqrt( gP' * P * gP ) = sqrt( g' * gP )
%   res_inf = ||g||_inf
%   step = || xnew - xold ||_P = || alpha * s ||_P 
%        = alpha * sqrt( s' * P * s )
%   step = || xnew - xold ||_inf
%
while (( (res > opts.g_tol) && (res_inf > opts.g_tol_inf) ) ...
    || ( (step > opts.x_tol) && (step_inf > opts.x_tol_inf) ) ...
    || ( f_old - f > opts.f_tol )) ...
    && (info.errno == 0)
  
  info.n_it = info.n_it + 1;

  % visualize current state
  if ~isempty(opts.visualize)
    opts.visualize(x);
  end
  
%% ******************* TRUST REGION SUBPROBLEM *********************
  
  % VERSION 1
  % NOTES:
  % minimize  q(s) = f + g^T s + 1/2 s^T P s  :  s^T P s <= DELTA^2
  % => (1-lambda) P s = - g
  % => s = - a * gP
  % q(a) = - a g^T gP + a^2/2  gP^T P gP = (-a+a^2/2) g^T gP
  % which is minimized for a = 1. So the minimizer subject to 
  % the constraint   a^2 g'*gP <= DELTA^2 is given by 
  a = min(1, DELTA/sqrt(g_dot_gP));
  s = - a * gP;
  % Note: g'*s = - a * g_dot_gP, s'*P*s = a^2 * g_dot_gP


%% ******************* UPDATES *********************

  % evaluate energy only
  [f_new, g_new] = ffun(x + s);
  info.n_f = info.n_f + 1;
  % check whether new state should be accepted
  rho = (f - f_new) / ((a-a^2/2) * g_dot_gP);
  if rho < opts.rho_tol    
    re = 'r';
    info.n_reject = info.n_reject + 1;
    % decrease trust region radius: DELTA_new = trial_step / 2.
    DELTA = a * sqrt(g_dot_gP) / 2;
    
   %     ACCEPT CASE
  else   
    re = 'a';
    % update states
    % g_old = g; gP_old = gP; x_old = x;
    g_dot_gP_old = g_dot_gP;
    f_old = f; 
    g = g_new; 
    f = f_new;
    x = x+s;
    % update preconditioner
    %     P = eval_P(Pfun, x, opts);
    %     gP = apply_P(P, g, opts);
    if opts.apply_P
      gP = Pfun(x, g);
    else
      % compute preconditioner
      P = eval_P(Pfun, x, opts);
      % precondition the first gradient
      gP = apply_P(P, g, opts);
    end
    info.n_p = info.n_p + 1;
    g_dot_gP = g'*gP;
    % update information
    res = sqrt(g_dot_gP);
    res_inf = norm(g, inf);
    step = sqrt( a^2 * g_dot_gP );  % s = - a * gP;
    step_inf = norm(s, inf);

    
    % update DELTA
    DELTA_old = DELTA;
    factor = min(0.5, (1 - DELTA_old / (2*sqrt(g_dot_gP_old))));
    DELTA = 2/factor * abs(f - f_old) / sqrt(g_dot_gP_old);
    DELTA = min(4*DELTA_old, max( DELTA_old / 4, DELTA ));
  end
  
  % display iteration information
  if opts.disp >= 1
    fprintf('%5d%s %5d %11.3e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n', ...
      info.n_it, re, info.n_f, f, res, res_inf, step, step_inf, ...
      DELTA, rho);
  end

  
  
%% ************* SAFEGUARDS *******************
  
  % bound on number of function evaluations
  if info.n_f > opts.maxn_f
    info.errno = 1;
    info.errmsg = 'n_f > maxn_f';
  end
  % lower bound on DELTA
  if DELTA < opts.DELTA_min
    info.errno = 2;
    info.errmsg = 'DELTA < DELTA_min';
  end


end



if opts.disp >= 1
  fprintf(['-----------------------------------------------------' ...
           '--------------------------------\n']);
  if info.errno <= 0
    fprintf(' popttr was succesful.\n');
  else
    fprintf([' popttr terminated with message ''', info.errmsg, '''\n']);
  end
  fprintf('   # Iterations                  : %d \n', info.n_it);
  fprintf('   # Function Evaluations        : %d \n', info.n_f);
  fprintf('   # Preconditioner Applications : %d \n', info.n_p);
  fprintf('   # Rejects                     : %d \n', info.n_reject);
  
  fprintf(['-----------------------------------------------------' ...
           '-----------------------\n']);
end



end



%% Auxiliary Routines for Preconditioning

% evaluate the preconditioner
% ~ is a placeholder for opts
function P = eval_P(Pfun, x, ~)
% case 1: Pfun is actually a function handle: evaluate P(x)
if isa(Pfun, 'function_handle')
  P = Pfun(x);
% case 2: Pfun is empty: set P = identity
elseif isempty(Pfun)
  P = 1;
% case 3: Pfun is a constant matrix - just return it
else
  P = Pfun;
end
end

% apply the preconditioner
function gP = apply_P(P, g, opts)
switch opts.Psolve
  case 'direct'
	gP = P \ g;
  
  case 'amg'
	hsl_mi20_startup;
	control = hsl_mi20_control;
	hsl_mi20_setup(P,control);  % inform = 
	gP = pcg(P,g,1e-8,2,'hsl_mi20_precondition');
	hsl_mi20_finalize;  % inform = 

  otherwise
	error('poptls: unknown value for opts.Psolve');
end
end






