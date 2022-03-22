%
% [a1, f1, g1, info] = wolfe(ffun, x, f, g, s, a1, scalf, opts, info)
%
% Strong Wolfe Linesearch for PDESCENT.
%
% Bracketing linesearch with cubic interpolation using function and
% gradient information at every step. Implementation following
% Nocedal & Wright, 2nd ed., Alg. 3.5 and 3.6.
% If the linesearch is unsuccesful and if the opts.switch_ls flag is set
% then it reverts to ls_armijo.
%

%
% List of Error messages:
%    0 : success
%  100 : 'alpha > alpha_max'
%  101 : 'existSW criterion failed'
%  102 : 'alpha < alpha_min'
%  103 : (free)
%  104 : 'insufficient progress in linesearch'
%  105 : '|a1-a0| < min_da'
%
function [a1, f1, g1, info] = wolfe( ...
  ffun, x, f, g, s, a1, scalf, opts, info)

% [Cubic linesearch with strong Wolfe conditions]
% computes alpha1, f1, g1 : x + alpha1 s satisfies strong Wolfe.
% based on Alg. 3.5 and 3.6 in (Nocedal & Wright, 2nd ed.) >> [NW]
% with some components taken from (More & Thuente, ACM, 1994) >> [MT]

% STEP 1: identify intercal with endpoints a0, a1,  containing an
%         admissible length. [NW, Alg. 3.5]

% store initial length
a_init = a1;
% directional derivative
d = s' * g;
% alpha_min and alpha_max are safeguards
amin = eps * scalf;
amax = opts.C1 * abs(f - opts.fmin) / abs(d);
% min_da is the smallest bracket size that is allowed
min_da = 100 * eps * scalf;
% suff_progress is the required progress after X iterations
% of the linesearch before it returns with an error message
suff_progress = max(1e-8 * scalf, opts.f_tol * 1e-2);
% the bracket is initialized with [a0, a1], a0 = 0, a1 is the
% user choice + the following safeguard
a0 = 0; f0 = f; g0 = g; d0 = d;
a1 = max(amin, min(a1, amax));

% evaluate \phi(\alpha_1) and \phi'(\alpha_1)
[f1, g1] = ffun(x + a1 * s);      
info.n_f = info.n_f + 1;
if (info.n_f > opts.maxn_f) || (f1 < opts.fmin), return; end
d1 = g1' * s;

% main loop of the bracketing phase
ls_i = 1;
while true
  if (f1 > f + opts.C1 * a1 * d) || ( (ls_i > 1) && (f1 >= f0) )
	% no bracketting needed, go to zoom stage with parameters (a0, a1)
    break;
  end
    
  if abs(d1) <= opts.C2 * abs(d)
    % a1 satisfies [A] (previous test) as well as [SW] (this test) so we 
	% don't need to do the "zoom" part. a1 alread contains the result!
    return;
  end
    
  if d1 >= 0
    % do zoom with reversed parameters (a1, a0)
    a2 = a1; f2 = f1; g2 = g1; d2 = d1;
    a1 = a0; f1 = f0; g1 = g0; d1 = d0;
    a0 = a2; f0 = f2; g0 = g2; d0 = d2;
    break;
  end
  
  % safeguard: if a1 == amax aleady, we need to stop now
  %            (this doesn't ever really happen though it seems)
  if (a1 == amax)
    if opts.switch_ls
      [a1, f1, g1, info] = ls_armijo( ...
        ffun, x, f, g, s, a_init, scalf, opts, info);
    else
      info.errno = 100;
      info.errmsg = 'alpha > alpha_max';
    end
    return;
  end
    
  % get new a1, initially stored in a2 as a0 becomes a1 now.
  a2 = min(a1 + 4 * (a1 - a0), amax);
  a0 = a1; f0 = f1; g0 = g1; d0 = d1;
  a1 = a2;
  
  % evaluate f, g, and \phi' at a1
  [f1, g1] = ffun(x + a1 * s);  
  info.n_f = info.n_f + 1;
  if (info.n_f > opts.maxn_f) || (f1 < opts.fmin), return; end
  d1 = g1' * s;
  
  ls_i = ls_i + 1;
end
  
% STEP 2: Zoom to find the admissible length
%         (Nocedal & Wright, 2nd ed., Alg. 3.6)
%

% store the previous 10 f values
nPrev = 10;
previous_f = inf * ones(nPrev, 1);
previous_f(1) = f1; previous_f(2) = f;

while true
  
  % debugging check: if neither of these three conditions
  % are satisfied then the algorithm has not worked!
  % (this doesn't ever seem to occur though!)
  existSW = false;
  if (f1 > f + opts.C1 * a1 * d) || (f1 >= f0) || (d1 >= 0)
    existSW = true;
  end
  if ~existSW
    if opts.switch_ls
      [a1, f1, g1, info] = ls_armijo( ...
        ffun, x, f, g, s, a_init, scalf, opts, info);
    else      
      info.errno = 101;
      info.errmsg = 'existSW criterion failed';
    end
    return;
  end
  
  % compute a safe-guardeed trial alpha, using cubic interpolation
  at = cubic_min(a0, f0, d0, a1, f1, d1);
  % safe-guard this step:
  bdry = abs(a1-a0) * 0.05;
  if isnan(at) || (at > max(a0,a1)-bdry) || (at < min(a0, a1)+bdry)
    at = 0.66 * a1 + 0.34 * a0;
  end
  ls_i = ls_i+1;
  
  % make sure things don't go crazy here
  if at < amin
    if opts.switch_ls
      [a1, f1, g1, info] = armijo( ...
        ffun, x, f, g, s, a_init, scalf, opts, info);
    else      
      info.errno = 102;
      info.errmsg = 'alpha < alpha_min';
    end
    return;
  end
  
  % evaluate f, g (TODO: should only evaluate f here!!)
  [ft, gt] = ffun(x+at * s);
  info.n_f = info.n_f + 1;
  if (info.n_f > opts.maxn_f) || (ft < opts.fmin), return; end
  dt = gt' * s;
  
  if opts.disp >= 2
    disp(['L.S. it. ', num2str(ls_i), ': alpha = ', num2str(at), ...
      ', [A] = ', num2str(ft - f - opts.C1 * at * d), ...
      ', [W] = ', num2str(abs(dt) - opts.C2 * abs(d)) ]);
  end
  
  % check whether suffient progress is being made in linesearch
  if ~(ft <= max(previous_f) - suff_progress)
    if opts.switch_ls
      [a1, f1, g1, info] = ls_armijo( ...
        ffun, x, f, g, s, a_init, scalf, opts, info);
    else      
      info.errno = 104;
      info.errmsg = 'insufficient progress in linesearch';
    end
    return;
  end
  % update the list of last 10 f.
  previous_f = [ft; previous_f(1:nPrev)];
  
  
  % First test for convergence (Armijo)
  if (ft > f + opts.C1 * at * d) || (ft >= f0)
	% required update if Armijo condition fails
    a1 = at; f1 = ft; g1 = gt; d1 = dt;
  else
    % Second test for convergence (Wolfe)
    if abs(dt) <= opts.C2 * abs(d)
      a1 = at; f1 = ft; g1 = gt; % d1 = dt;
      break;
    end
    % no convergence yet >> continue testing
    if dt * (a1 - a0) >= 0
      a1 = a0; f1 = f0; g1 = g0; d1 = d0;
      a0 = at; f0 = ft; g0 = gt; d0 = dt;
      if opts.disp >= 2
        disp(['Update: a0 = ', num2str(a0)])
      end
    else
      a0 = at; f0 = ft; g0 = gt; d0 = dt;
      if opts.disp >= 2
        disp(['Update: a0 = ', num2str(a0)])
      end
    end
  end
  
  if abs(a1-a0) < min_da
    if opts.switch_ls
      [a1, f1, g1, info] = ls_armijo( ...
        ffun, x, f, g, s, a_init, scalf, opts, info);
    else      
      info.errno = 105;
      info.errmsg = '|a1-a0| < min_da';
    end
    % again: switch to backtracking?
    return;
  end
end

end


%
% Auxiliary function: compute minimizer of a cubic polynomial
% will return NaN if calculation is unreliable.
%
function ac = cubic_min(a0, f0, d0, a1, f1, d1)
% (h1, h2 are d1, d2 in (3.59))
h1 = d0 + d1 - 3 * (f0 - f1) / (a0 - a1);
if h1^2 - d0*d1 <= 1e-10 * abs(a1-a0)
  ac = NaN;
  return;
else
  h2 = sign(a1 - a0) * sqrt(h1^2 - d0 * d1);
  if abs(d1-d0+2*h2) <= 1e-8 * abs(d1+h2-h1)
    ac = NaN;
  else
    ac = a1 - (a1-a0) * ((d1+h2-h1) / (d1-d0+2*h2));
  end
end
end


% %
% % Auxiliari function: compute minimizer of quadratic
% % returns NaN if calculation is unreliable
% %
% function ac = quad_min_frc(a0, f0, d0, a1, f1, d1)
% 
% end

