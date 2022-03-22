% [x, info] = popttr1(x, ffun, P, H, options)
%
% Preconditioned trust region method, basic Steihaug method
% with global rank-1 updates of the hessian
%
% Input:
%   x0       : initial guess
%   ffun     : objective function; [f, g] = ffun(x)
%   P        : preconditioner class
%   H        : approximate hessian class
%   options  : structure (see popttr_options for details)
%
% Output:
%   x    : final iterate, empty array if psd was unsuccesful
%   info : structure with the following content
%            .n_it : number of iterations
%            .n_f  : number of function evaluations
%            .n_p  : number of preconditioner applications
%            .errmsg
%            .n_reject


% INTERFACE FOR P

% INTERFACE FOR H


%          if opts.hist = true, the following history is stored in info:
%            .h_res, .h_res_inf, .h_step, .h_step_inf
%            .h_f, .h_n_f, .h_n_it


function [x, info] = popttr1(x, ffun, P, H, opts)


%% *********************  Initialization  *************************

% initialize info structure
info.errno = 0; 
info.errmsg = 'success';
info.n_it = 0; 
info.n_f = 0; 
info.n_p = 0;
info.n_h = 0;
info.n_reject = 0;
% first function evaluation
[f, g] = ffun(x);
info.n_f = info.n_f + 1;
% compute preconditioner
P = P.init(P);
P = P.eval(P, x);
% Hessian initialisation
H = H.init(H, P);
% precondition the first gradient
gP = P.apply(P, g);
info.n_p = info.n_p + 1;
% auxiliary scalar: g'*gP = g'*P^{-1}*g = gP'*P*gP
g_dot_gP = g'*gP;
% initialize trust region radius
DELTA = min(1, sqrt(abs(gP'*g)));

%% ********************* Main Loop ******************
res = sqrt(g_dot_gP); res_inf = norm(g, inf);
step = 0; step_inf = 0;
f_old = f + opts.f_tol + 1;
if opts.disp >= 1
    fprintf(['-----------------------------------------------------' ...
             '--------------------------------\n']);
    fprintf('%5s %11s %9s %9s %9s %9s %9s %9s\n', ...
      'n_it', 'f', 'Res(P)', 'Res(inf)', 'Step(P)', 'Step(Inf)', ...
      'DELTA', 'rho');
    fprintf(['-----------------------------------------------------' ...
             '--------------------------------\n']);
    fprintf('%5d %11.3e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n', ...
      info.n_it, f, res, res_inf, step, step_inf, DELTA, 0);
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
  
  % minimize  q(s) = f + g^T s + 1/2 s^T H s  :  s^T P s <= DELTA^2
  % using the Steihaug method
  s = steihaug(g, H, P, DELTA);
  

  % get sxg, sxHxs, step = sqrt(s'*P*s)

%% ******************* UPDATES *********************

  % evaluate energy only
  [f_new, g_new] = ffun(x + s);
  info.n_f = info.n_f + 1;
  % check whether new state should be accepted
  rho = (f - f_new) / (- sxg - 0.5 * sxHxs);
  if rho < 0.01    
    re = 'r';
    info.n_reject = info.n_reject + 1;
    % decrease trust region radius: DELTA_new = trial_step / 2.
    DELTA = step / 2;
    
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
    P = P.eval(P, x);
    gP = P.apply(P, g);
    info.n_p = info.n_p + 1;
    g_dot_gP = g'*gP;
    % update approximate hessian
    H = H.update(H, P, x, g, gP);
    % update information
    res = sqrt(g_dot_gP);
    res_inf = norm(g, inf);
    step_inf = norm(s, inf);
    % step was already updated during steihaug

    % update DELTA  TODO TODO TODO TODO TODO TODO
    DELTA_old = DELTA;
    DELTA = (1 - DELTA_old / (2*sqrt(g_dot_gP_old)))^(-1) * ...
      2 * abs(f - f_old) / sqrt(g_dot_gP_old);
    DELTA = min(4*DELTA_old, max( DELTA_old / 4, DELTA ));
  end
  
  % display iteration information
  if opts.disp >= 1
    fprintf('%5d%s %11.3e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n', ...
      info.n_it, re, f, res, res_inf, step, step_inf, ...
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
    fprintf([' pdescent terminated with message ''', info.errmsg, '''\n']);
  end
  fprintf('   # Iterations                  : %d \n', info.n_it);
  fprintf('   # Function Evaluations        : %d \n', info.n_f);
  fprintf('   # Preconditioner Applications : %d \n', info.n_p);
  fprintf('   # Rejects                     : %d \n', info.n_reject);
  
  fprintf(['-----------------------------------------------------' ...
           '-----------------------\n']);
end



end




%% STEIHAUG METHOD for TR SUBPROBLEM -- SIMPLIFIED VERSION
% g = gradient, v = P^{-1} g
%
% return: s, nH, nP
%
function [s, nrm_s] = steihaug(g, v, P, H)

% initialize vectors
s = zeros(size(g));
p = - v;
Hp = H.apply(H, p);

% initialize inner products
g_dot_v = g' * v;
sPs = 0; 
sPp = 0; 
pPp = g' * v;
pHp = p' * Hp;
sHp = 0;
sg = 0;

while true
  % check for negative curvature ...
  if pHp <= 0
    % and if this occurs, just move to the boundary of the trust region 
    % and then stop:
    % find tau s.t. ||s + tau * p||_P = DELTA
    %  tau^2 pPp + 2 tau sPp + sPs - DELTA^2 = 0
    tau = (- sPp + sqrt(sPp^2 + DELTA^2 - sPs)) / pPp;
    s = s + tau * p;
    nrm_s = DELTA;
    return;
  end
  
  % have positive curvature direction => try standard CG step
  alpha = g_dot_v / (pHp+eps*g_dot_v);
  % check whether this step would take us outside ... 
  if (sPs + 2*alpha * sPp + alpha^2 * pPp > DELTA^2)
    % ... and if so cut it at the boundary, and then return
    % (same as above)
    tau = (- sPp + sqrt(sPp^2 + DELTA^2 - sPs)) / pPp;
    s = s + tau * p;
    nrm_s = DELTA;
    return;
  end
  
  % standard CG updates for directions
  s = s + alpha * p;
  g = g + alpha * Hp;
  v = P.apply(P, g);
  % next line does the following:
  %       beta = gv_{new} / gv_{old} + compute gv_{new}
  beta = g_dot_v; g_dot_v = g' * v; beta = g_dot_v / beta;
  p = - v + beta * p;
  
  % update Hp
  Hp = H.apply(H, p);
  
  % update the inner products
  % first group
  sHs = sHs + 2 * alpha * sHp + alpha^2 * pHp;
  sg = sg + alpha * sHp;
  pHp = p' * hP;
  sHp = sHp + alpha * pHp;
  g_dot_v = g' * v;
  % second group
  sPs = sPs + 2 * alpha * sPp + alpha^2 * pPp;
  sPp = - sg + beta * sPp + alpha * beta * pPp;
  pPp = g_dot_v + beta^2 * pPp;
end

end



% 
% %% Auxiliary Routines for Preconditioning
% 
% % evaluate the preconditioner
% % ~ is a placeholder for opts
% function P = eval_P(Pfun, x, ~)
% % case 1: Pfun is actually a function handle: evaluate P(x)
% if isa(Pfun, 'function_handle')
%   P = Pfun(x);
% % case 2: Pfun is empty: set P = identity
% elseif isempty(Pfun)
%   P = 1;
% % case 3: Pfun is a constant matrix - just return it
% else
%   P = Pfun;
% end
% end
% 
% % apply the preconditioner
% function gP = apply_P(P, g, opts)
% switch opts.Psolve
%   case 'direct'
% 	gP = P \ g;
%   
%   case 'amg'
% 	hsl_mi20_startup;
% 	control = hsl_mi20_control;
% 	hsl_mi20_setup(P,control);  % inform = 
% 	gP = pcg(P,g,1e-8,2,'hsl_mi20_precondition');
% 	hsl_mi20_finalize;  % inform = 
% 
%   otherwise
% 	error('poptls: unknown value for opts.Psolve');
% end
% end
% 
% 
% 



