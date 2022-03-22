%
% [a1, f1, g1, info] = armijo(ffun, x, f, g, s, a1, scalf, opts, info)
%
% Armijo Linesearch for PDESCENT
%
% Backtracking algorithm which uses a cubic interpolant at each step,
% with safeguard mechanisms that ensure sufficient progress at each
% step.
%

%
%   202 : 'alpha < alpha_min'
%   204 : 'insufficient progress in linesearch'
%
function [a1, f1, g1, info] = armijo( ...
  ffun, x, f, g, s, a1, scalf, opts, info)

% directional derivative
d = s' * g;
% alpha_min and alpha_max are safeguards
amin = min(eps * 10 * scalf, opts.f_tol * 1e-2);
amax = opts.C1 * abs(f - opts.fmin) / abs(d);
% suff_progress is the required progress after X iterations
% of the linesearch before it returns with an error message
suff_progress = max(1e-8 * scalf, opts.f_tol * 1e-2);
% a1 is the user choice + the following safeguard
a1 = max(amin, min(a1, amax));

% evaluate \phi(\alpha_1) and \phi'(\alpha_1)
[f1, g1] = ffun(x + a1 * s);      
info.n_f = info.n_f + 1;
if info.n_f > opts.maxn_f, return; end
d1 = g1' * s;

% main loop of the bracketing phase
ls_i = 1;
% store the previous 10 f values
nPrev = 10;
previous_f = inf * ones(nPrev, 1);
previous_f(1) = f1; previous_f(2) = f;

while f1 > f + opts.C1 * a1 * d
  % compute a safe-guardeed trial alpha, using cubic interpolation
  at = cubic_min(0, f, d, a1, f1, d1);
  % safe-guard this step against nan and against too little 
  % or too much progress
  if isnan(at) || isinf(at)
    at = 0.5 * a1;
  end
  at = max( 0.1 * a1, min( at, 0.8 * a1 ) );
  ls_i = ls_i+1;
  
  % make sure things don't go crazy here
  if at < amin
    info.errno = 202;
    info.errmsg = 'alpha < alpha_min';
    return;
  end
  
  % evaluate f, g (TODO: should only evaluate f here!!)
  a1 = at;
  [f1, g1] = ffun(x+a1 * s);
  info.n_f = info.n_f + 1;
  if info.n_f > opts.maxn_f, return; end
  d1 = g1' * s;
  
  % check whether sufficient progress is being made
  if ~(f1 <= max(previous_f) - suff_progress)
    info.errno = 204;
    info.errmsg = 'insufficient progress in linesearch';
    return;
  end
  % update the list of last 10 f.
  previous_f = [f1; previous_f(1:nPrev)];
  
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
