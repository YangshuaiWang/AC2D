
% [x, info] = poptls(ffun, Pfun, x0, options)
% Preconditioned linesearch optimization algorithms, including
% steepest descent, nonlinear conjugate gradients, LBFGS.
%
% Input:
%   ffun     : objective function; [f, g] = ffun(x)
%   Pfun     : preconditioner; P = Pfun(x), or constant matrix P, or []
%   x0       : initial guess
%   options  : structure (see psd_options)
%   data     : structure which may contain global or temporary 
%              information and is passed through all procedured
%              the entry data.opt is reserved and will contain
%              information for the optimization.
% Output:
%   x    : final iterate, empty array if psd was unsuccesful
%   info : structure with the following content
%            .n_it : number of iterations
%            .n_f  : number of function evaluations
%            .n_p  : number of preconditioner applications
%          if opts.hist = true, the following history is stored in info:
%            .h_res, .h_res_inf, .h_step, .h_step_inf
%            .h_f, .h_n_f, .h_n_it


function [x, info] = poptls(ffun, Pfun, x, opts)

%% *********************  Initialization  *************************

if strcmp(opts.Psolve, 'amg')
  	addpath ~/Library/mTools/hsl_mi20
end

% initialize info structure
info.errno = 0;
info.errmsg = 'success';
info.n_it = 0; 
info.n_f = 0; 
info.n_p = 0;
info.n_restart = 0;
% first function evaluation
[f, g] = ffun(x);
info.n_f = info.n_f + 1;
% determine preconditioner (dynamic, none, static)
if isa(Pfun, 'function_handle')
  P = Pfun(x);
elseif isempty(Pfun)
  P = 1;
else
  P = Pfun;
end
% precondition the first gradient
% gP = P \ g;
gP = apply_P(P, g, opts);
info.n_p = info.n_p + 1;
% start the scaling of the problem
if isempty(opts.scalf)
  scalf = max(abs(f), sqrt(abs(g'*gP)));
  dyn_scalf = true;
else
  scalf = opts.scalf;
  dyn_scalf = false;
end
% the following choice of s_old guarantees that the first step
% taken is always a (preconditioned) steepest descent step.
s_old = 0 * g;  
g_old = g; gP_old = gP; 
% force steepest descent step
force_sd = false;
force_armijo = false;
% restart_done is to tell the output whether the current step
% had to perform a restart of CG method
restart_done = false;


%% ********************* Main Loop ******************
res = sqrt(g' * gP);
res_inf = norm(g, inf);
step = 0;
step_inf = 0;
f_old = f + opts.f_tol + 1;
if opts.disp >= 1
    fprintf(['-----------------------------------------------------' ...
             '-----------------------\n']);
    fprintf('%5s %6s %11s %9s %9s %9s %9s %9s \n', ...
      'n_it', 'n_f', 'f', 'Res(P)', 'Res(inf)', 'Step(P)', 'Step(Inf)', ...
      'alpha');
    fprintf(['-----------------------------------------------------' ...
             '-----------------------\n']);
    fprintf('%5d %6d %11.3e %9.2e %9.2e %9.2e %9.2e %9.2e \n', ...
      info.n_it, info.n_f, f, res, res_inf, step, step_inf, 0.0);
end

% Notes on the termination criterion:
%   res = sqrt( gP' * P * gP ) = sqrt( g' * gP )
%   res_inf = ||g||_inf
%   step = || xnew - xold ||_P = || alpha * s ||_P 
%        = alpha * sqrt( s' * P * s )
%   step = || xnew - xold ||_inf
%
while ( (res > opts.g_tol) && (res_inf > opts.g_tol_inf) ) ...
    || ( (step > opts.x_tol) && (step_inf > opts.x_tol_inf) ) ...
    || ( f_old - f > opts.f_tol )
  
  % save iteration history  
  if opts.hist
    info.h_n_it(info.n_it+1) = info.n_it;
    info.h_n_f(info.n_it+1) = info.n_f;
    info.h_n_p(info.n_it+1) = info.n_p;
    info.h_res(info.n_it+1) = res;
    info.h_step(info.n_it+1) = step;
    info.h_res_inf(info.n_it+1) = res_inf;
    info.h_step_inf(info.n_it+1) = step_inf;
    info.h_f(info.n_it+1) = f;
  end
  
  % visualize current state
  if ~isempty(opts.visualize)
    opts.visualize(x);
  end
  
  %% ******************* Choice of Search Direction *****************
  % Compute the beta-parameter for the choice of search direction
  % Nocedal & Wright, 2nd ed., Sec. 5.2
  % TODO: HS method, (5.49) and (5.50) methods
  if force_sd
    beta = 0; %#ok
    s = - gP;
    d = g' * s;
    force_sd = false;
    
  elseif strcmp(opts.direction, 'lbfgs')
    if info.n_it > 0
      % update
      lbfgs_data = lbfgs_update(lbfgs_data, x - x_old, g - g_old);
      s = lbfgs_apply(lbfgs_data, -g);
      d = s' * g;
    else
      % initialize
      lbfgs_data = lbfgs_init(x, P, opts);
      s = -g;
      d = g' * s;
    end
    
    % if s is not a descent direction, restart the LBFGS
    if (d >= - eps * 10 * scalf)
      lbfgs_data = lbfgs_init(x, P, opts);
      s = -g;
      d = g' * s;
      info.n_restart = info.n_restart + 1;
      restart_done = true;
    end
      
    
  else
    switch opts.direction
      case 'sd'    % Steepest Descent
        beta = 0;
      case 'fr'    % Fletcher Reeves
        beta = (gP' * g) / (gP_old' * g_old);
      case 'pr'    % Original Polak-Ribiere (no convergence guarantee!)
        beta = (gP'*(g - g_old)) / (gP_old' * g_old);
      case 'pr+'   % Polak-Ribiere, safe-guarded by SD
        beta_pr = (gP'*(g - g_old)) / (gP_old' * g_old);
        beta = max(0, beta_pr);
      case 'fr-pr' % Polak-Ribiere, safe-guarded by FR
        beta_fr = (gP' * g) / (gP_old' * g_old);
        beta_pr = (gP'*(g - g_old)) / (gP_old' * g_old);
        beta = max( -beta_fr, min( beta_fr, beta_pr ) );
      otherwise    % Unkown method
        error(['plinesearch: Unknown search direction: ', opts.direction]);
    end
    
    % compute the (conjugate) search direction
    s = - gP + beta * s_old;
    % Safeguard: Restart with steepest descent if
    %   - the directional derivative is (numerically) positive,
    %   - or successive gradients are far from orthogonal
    d = g' * s;
    if (beta ~= 0) && ...
        ( (d >= - eps * 10 * scalf) ...
        || (abs((gP' * g_old) / (gP' * g)) > 0.8) )
      beta = 0; %#ok
      s = -gP;
      d = g' * s;
      info.n_restart = info.n_restart + 1;
      restart_done = true;
    end
  end
  
  % if the directional derivative too large then stop
  if (d >= - eps * 100 * scalf)
    info.errno = 20;
    info.errmsg = 'directional derivative >= 0 (numerically)';
    break;
  end

  
  %% ********* Selection of initial steplength for linesearch ********
  % the following choices will be modified by some safeguards below
  if info.n_it > 0
    switch opts.ls_init
	  case -1
		% for pure Newton: force initial trial step of 1
		% (more for test purposes than anything else)
		a1 = 1;
		
      case 0
        % for Newton, quasi-Newton, etc.
        a1 = min(1, 3 * abs(f - f_old) / abs(g_old' * s_old));

      case 1
        % extrapolated from the previous iteration so that [A] would
        % holds with C1 = 1/2, which is somehow optimal and would also
        % give superlinear convergence
        a1 = 2 * abs(f - f_old) / abs(g_old' * s_old);
      case 2
        % similar to 1
        a1 = 2 * abs(f - f_old) / abs(g' * s);
      case 3
        a1 = alpha_old * abs( (g_old' * s_old) / (g' * s) );
      otherwise
        error('plinesearch: unknown choice of opts.ls_init');
    end
    a1 = min(4*alpha_old, max( alpha_old / 4, a1 ));
  else 
    % on the first iteration, make a different choice!
    a1 = min([opts.alpha_init, 1, 1 / sqrt(abs(gP'*g))]);
  end
      
 
  %% ******************** Linesearch ************************
  % after the linesearch, we  will have obtained a1 satisfying the 
  % Armijo and/or strong Wolfe condition, and  
  % [f1, g1] = ffun(x+a1*s) will be evaluated.
  if force_armijo || strcmp(opts.linesearch, 'armijo');
    [a1, f1, g1, info] = armijo( ...
      ffun, x, f, g, s, a1, scalf, opts, info);
    force_armijo = false;    
  elseif strcmp(opts.linesearch, 'wolfe')
    [a1, f1, g1, info] = wolfe( ...
      ffun, x, f, g, s, a1, scalf, opts, info);
    % if switching to armijo is allowed and there are certain
    % types of errors, restart with a steepest descent step
    if opts.force_armijo && (info.errno >= 100)
      force_sd = true;
      force_armijo = true;
      info.errno = 0;
      info.errmsg = 'success';
      info.n_restart = info.n_restart + 1;
      restart_done = true;
      continue;
    end
  else
      error('plinesearch: Unknown linesearch option.');
  end
  % check whether an error occured during linesearch
  if info.errno > 0
    break;
  end
  
  %% ***************** Do all the updates ************************
  % store the old gradient, preconditioned gradient, direction
  g_old = g; gP_old = gP; s_old = s; f_old = f; x_old = x;
  % update energy, gradient, preconditioned gradient
  if isa(Pfun, 'function_handle')
    P = Pfun(x);
  end
  f = f1; g = g1; 
  % gP = P \ g;          
  gP = apply_P(P, g, opts);
  info.n_p = info.n_p + 1;
  % new iterate
  x = x + a1 * s;                      
  info.n_it = info.n_it + 1;
  % the residual is computed at the new position
  res = sqrt(g' * gP); 
  res_inf = norm(g, inf);
  % the steplength taken in this iteration
  step = a1 * sqrt(s' * P * s);
  step_inf = a1 * norm(s, inf);
  % the old alpha_value
  alpha_old = a1;
  % rescale the problem
  if dyn_scalf
    scalf = max([scalf, abs(f), sqrt(abs(g'*gP))]);
  end
  
  if opts.disp >= 1
    if restart_done
      re = 'r';
    else
      re = ' ';
    end
    fprintf('%5d%s %5d %11.3e %9.2e %9.2e %9.2e %9.2e %9.2e \n', ...
      info.n_it, re, info.n_f, f, res, res_inf, step, step_inf, a1);
    restart_done = false;
  end
  
  % check for erroneous termination:
  info = check_termination(f, info, opts);
  if info.errno > 0, break; end
end

if opts.disp >= 1
  fprintf(['-----------------------------------------------------' ...
           '-----------------------\n']);
  if info.errno <= 0
    fprintf(' pdescent was succesful.\n');
  else
    fprintf([' pdescent terminated with message ''', info.errmsg, '''\n']);
  end
  fprintf('   # Iterations                  : %d \n', info.n_it);
  fprintf('   # Function Evaluations        : %d \n', info.n_f);
  fprintf('   # Preconditioner Applications : %d \n', info.n_p);
  fprintf('   # Restarts                    : %d \n', info.n_restart);
  
  fprintf(['-----------------------------------------------------' ...
           '-----------------------\n']);
end

end



%
% Auxiliary function which outputs an error if the 
% procedure doesn't terminate.
%
function info = check_termination(f, info, opts)
if f < opts.fmin
  info.errno = 1;
  info.errmsg = 'failure because f < fmin';
  return;
end
if info.n_f > opts.maxn_f
  info.errno = 2;
  info.errmsg = 'failure because n_f > max_f';
  return;
end
if info.n_it > opts.maxn_it
  info.errno = 3;
  info.errmsg = 'failure because n_it > maxn_it';
  return;
end
end



function q = lbfgs_apply(lbfgs_data, q)
%
% see [NW, Alg. 7.4]
%
a = zeros(1, lbfgs_data.Mact);
for m = 1:lbfgs_data.Mact
  a(m) = lbfgs_data.RHO(m) * lbfgs_data.S(:,m)' * q;
  q = q - a(m) * lbfgs_data.Y(:, m);
end
gamma = (lbfgs_data.S(:, 1)' * lbfgs_data.Y(:, 1)) ...
  / (lbfgs_data.Y(:, 1)' * lbfgs_data.Y(:, 1));
q = gamma * (lbfgs_data.B0 \ q);
% q = gamma * q;
for m = lbfgs_data.Mact:(-1):1
  beta = lbfgs_data.RHO(m) * lbfgs_data.Y(:, m)' * q;
  q = q + lbfgs_data.S(:, m) * (a(m) - beta);
end
end

function lbfgs_data = lbfgs_update(lbfgs_data, s_new, y_new)
Mmax = lbfgs_data.Mmax;
% shift the previous vectors to the back
lbfgs_data.S(:, 2:Mmax) = lbfgs_data.S(:, 1:(Mmax-1));
lbfgs_data.Y(:, 2:Mmax) = lbfgs_data.Y(:, 1:(Mmax-1));
lbfgs_data.RHO(2:Mmax) = lbfgs_data.RHO(1:(Mmax-1));
% add in the new vectors
lbfgs_data.S(:, 1) = s_new;
lbfgs_data.Y(:, 1) = y_new;
lbfgs_data.RHO(1) = 1 / (s_new' * y_new);
% update M-data
lbfgs_data.Mact = min(lbfgs_data.Mact+1, lbfgs_data.Mmax);
end

function [lbfgs_data] = lbfgs_init(x, P, opts)
lbfgs_data.Mmax = opts.lbfgs_Mmax;
lbfgs_data.Mact = 0;
lbfgs_data.S = zeros(numel(x), lbfgs_data.Mmax);
lbfgs_data.Y = zeros(numel(x), lbfgs_data.Mmax);
lbfgs_data.RHO = zeros(1, lbfgs_data.Mmax);
lbfgs_data.B0 = P;
end



function gP = apply_P(P, g, opts)
switch opts.Psolve
  case 'direct'
	gP = P \ g;
  
  case 'amg'
	hsl_mi20_startup;
	control = hsl_mi20_control;
	inform = hsl_mi20_setup(P,control);
	gP = pcg(P,g,1e-8,2,'hsl_mi20_precondition');
	inform = hsl_mi20_finalize;

  otherwise
	error('poptls: unknown value for opts.Psolve');
end

end




