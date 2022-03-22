
% opts = poptls_opts(direction, Pquality)
%   Generates default options for the preconditioned descent routine
%   PLINESEARCH. It accepts two optional parameters:
%
% direction : as below: 'sd', 'fr', 'pr', 'pr+', 'fr-pr', or 'lbfgs'
% Pquality  : quality of the preconditioner. possible choices are
%             'none' : problem may be poorly conditioned, and no
%                      preconditioner is used
%            'rough' : the topology induced by the preconditioner is
%                      is roughly correct, but little or no Hessian 
%                      information is included
%             'good' : between 'rough' and 'Newton'
%           'newton' : Hessian, or very close.            
% the choice of Pquality affects how the linesearch parameters are
% chosen.
%
% Primary method parameters
%   .direction : Which steepest descent or conjugate gradient method:
%            'sd' : steepest descent
%            'fr' : Fletcher-Reeves
%            'pr' : Polak-Ribiere
%            'pr+' : safe-guarded modified PR (take max. with 0)
%            'fr-pr' : max( -fr, min(fr, pr) )
%            'lbfgs' : limited memory bfgs (no preconditioning!)
%            *** TODO: Hestenes/Stiefel (5.46), (5.49), (5.50) ***
%   .linesearch : 'wolfe' or 'armijo'
%   .ls_init : how to initialize the linesearch (subject to safeguards):
%              0 : alpha0 = min(1, 2.2 (f-fold) / (g_old * s_old))
%             {1}: alpha0 ~ 2 (f-fold) / (g_old * s_old)
%              2 : alpha0 ~ 2 (f-fold) / (g * s)
%              3 : alpha0 ~ alpha_old * (g_old * s_old) / (g * s);
%   .Psolve : 'direct' (backslash), 'amg'
%
% Termination tolerance, succesful termination requires three criteria:
%   .g_tol, .g_tol_inf : for succesful termination need 
%                          either || P^{-1}g ||_P < g_tol
%                          or || g ||_inf < g_tol_inf
%   .x_tol, .x_tol_inf : for succesful termination need 
%                          either || xnew - x ||_P < x_tol
%                          or || xnew - x ||_inf < x_tol_inf
%   .f_tol             : for succesful termination need 
%                          |f_old - f_new| < f_tol
%
% Unsuccesful termination
%   .maxn_nit : maximum number of iterations
%   .maxn_f   : maximum number of fcn evaluations
%   .fmin     : if f goes below fmin then the algorithm terminates
%
% Output parameters
%   .hist : {false} If true, then the info field contains information 
%           on the history of the iteration.
%   .disp : level of detail for console output
%           {0}: after termination, 1: after each iteration, 2: debugging
%   .visualize : function handle visualize(x), that is called at the
%                beginning of every iteration
%
% Linesearch Parameters
%   .C1 : linesearch parameter for Armijo condition
%   .C2 : linesearch parameter for the strong Wolfe condition
%   .switch_ls : if true, then an unsuccesful wolfe ls can switch to
%                armijo linesearch (default is {false})
%   .force_armijo : forces one armijo linesearch after every unsuccesful
%                   iteration, when the direction is reset to the steepest
%                   descent direction. (default is {true})
%
% LBFGS Parameters
%   .lbfgs_store : number of past iterates to store in the LBFGS method
%
% Scaling Parameters
%   .alpha_init : first trial steplength
%   .scalf : a constant that determines a natural scaling for the
%            problem. This is normally set automatically, but can 
%            be enforced here. (default {[]})
%

function opts = poptls_opts(direction, Pquality)
% Default Parameters
if nargin < 1, direction = 'pr+'; end
if nargin < 1, Pquality = 'none'; end
% these will be set at the end. The following are default parameters
% that will in fact be overwritten below.

%% The three main choices in the design of the method are the following:
% Which descent direction to take?
opts.direction = 'pr+';
% Which line search {'armijo', 'wolfe'}
opts.linesearch = 'wolfe';
% How to choose the initial alpha?
opts.ls_init = 1;
% how to invert the preconditioner
opts.Psolve = 'direct';

%% Termination criteria
% criteria for succesful termination
opts.g_tol = 1e-4;
opts.g_tol_inf = 1e-6;
opts.x_tol = 1e-5;
opts.x_tol_inf = 1e-7;
opts.f_tol = 1e-6;
% Criteria for early termination
opts.maxn_it = 1e3;
opts.maxn_f = 1e4;
opts.fmin = -1e300;

%% Linesearch Parameters
% Parameter for Armijo condition
opts.C1 = 1e-4;
% Parameter for the Wolfe condition
opts.C2 = 0.9;
% if linesearch is wolfe, this flag allows to switch to armijo if
% a disaster occurs (switch from within ls_wolfe)
opts.switch_ls = false;
% similarly, the following flag allows a restart with steepest
% descent AND armijo if a Wolfe linesearch has failed.
opts.force_armijo = true;

%% output parameters
% if true, then details about the iteration history are returned
opts.hist = false;
% display level: 0: no display, 1: after each iteration, 2: debugging
opts.disp = 0;
% if this is set to a function handle then this function is called
% after every iteration with parameters (x, varargin{:})
opts.visualize = [];

%% Some scaling constants
% initial step length. The actual initial steplength is a combination 
% of different things
opts.alpha_init = 1;
% a scaling factor which is dynamically updated if scalf = []
opts.scalf = [];

%% Options for the LBFGS method
opts.lbfgs_Mmax = 100;


%% Set default parameters for direction / Pquality
% TODO: ensure direction string is admissible
opts.direction = direction;

if strcmp(direction, 'lbfgs')
  opts.C2 = 0.9;
  opts.C1 = 1e-4;
  opts.ls_init = 0;
  opts.maxn_it = 300;
  opts.maxn_f = 1000;
else
  switch(Pquality)
    case 'none'
      opts.C2 = 0.2;
      opts.C1 = 1e-4;
      opts.ls_init = 3;
      opts.maxn_it = 1e3;
      opts.maxn_f = 1e4;
      
    case 'rough'
      opts.C2 = 0.5;
      opts.C1 = 1e-4;
      opts.ls_init = 1;
      opts.maxn_it = 300;
      opts.maxn_f = 1000;
      
    case 'good'
      opts.C2 = 0.7;
      opts.C1 = 1e-4;
      opts.ls_init = 0;
      opts.maxn_it = 300;
      opts.maxn_f = 1000;
      
    case 'newton'
      opts.C2 = 0.9;
      opts.C1 = 1e-4;
      opts.ls_init = 0;
      opts.maxn_it = 300;
      opts.maxn_f = 1000;
      
    otherwise
      error('plinesearch_options: illegal Pquality parameter');
  end
end

%% Trust region options
% if set to true, then Pfun _applies_ the preconditioner
% instead of generating it
opts.apply_P = false;
opts.DELTA_min = 1e-10;
opts.rho_tol = 1e-2;

end
