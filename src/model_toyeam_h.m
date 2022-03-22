
%% function model = model_toyeam_h(a, b, c, rH)
% generates a model structure for a toy EAM model
%    V(r) = sum_j J(|r_j|) + F( \sum_j \rho(|r_j|) )
%    where   * J(s) = exp(-2a(s-1)) - 2*exp(-a(s-1))
%            * \rho(s) = exp(-bs)
%            * F(t) = c((t - t0)^2 + (t - t0)^4), where t0 = 6 \rho(0.9).
% for this model the free variable is the deformation field.
%
% Input
%       a : stiffness/decay of pair potential
%       b : decay of electron density
%       c : stiffness of embedding energy
%      rH : cut-off in hopping distance
%
% Output
%   model : structure containing the following fields
%      * model.a, model.b, model.c   (input parameters)
%      * model.rCutH : 1      (just nearest neighbour model)
%      * model.id = 'toyeam:tri2d'
%      * model.Vfun : site potential, [V, dV] = Vfun(model, r)
%      * model.Wfun : CB potential [W, dW] = Wfun(model, F)
%      * model.rDim = 2, model.dDim = 2  (dimension of range and domain)
%      * model.A : transformation matrix for triangular lattice
%                  from Z^2 to A Z^2    (probably not required)
%      * model.rCB : set of basis vectors for computing W
%


% TODO



function model = model_toyeam_h(a, b, c, rH)

if nargin == 0
  test_model_toyeam();
  return;
end


model.a = a;
model.b = b;
model.c = c;
model.rCutH = rH;
model.rCut = [];
model.id = 'tri2d_toyeam_h';
model.Vfun = @toyeam_site_potential;
model.Wfun = @Wcb;
% the basic ground state should be the triangular lattice
model.A = [ [1;0] [cos(pi/3); sin(pi/3)] ];
model.rDim = 2;
model.dDim = 2;
model.bc = '';

% get vectors for computing the CB potential
model.rCB = get_nhd_hop(model.A, rH);

end


%% function [V, dV] = morse_site_potential(model, r)
% Input
%   model : structure, must contain model.a, model.r0, model.rCut
%   r     : d x n double array, d = dimension, n = number of neighbours
% Output
%   V  : double, site energy
%   dV : d x n double array, forces on the n neighbours
%
function [V, dV] = toyeam_site_potential(model, r)
% parameters
d = model.dDim;
% compute bond lengths
s = sqrt(sum(r.^2, 1));

% 1. evaluate the pair potential part  (nearest-neighbour only!)
if nargout == 1
  % energy only
  J = exp(-2*model.a*(s-1)) - 2 * exp(-model.a*(s-1));
  % J = MLJpotential(s, model.a, 1, []);
  V = 0.5 * sum(J);
else  
  J = exp(-2*model.a*(s-1)) - 2 * exp(-model.a*(s-1));
  dJ = -2*model.a * ( exp(-2*model.a*(s-1)) - exp(-model.a*(s-1)) );
  % [J, dJ] = MLJpotential(s, model.a, 1, []);
  % energy
  V = 0.5 * sum(J);
  % forces
  dV = 0.5 * (ones(d,1) * (dJ./s)) .* r;
end

% 2. evaluate the EAM part
rho = exp(- model.b * s);
t = sum(rho);
t0 = 6 * exp(-model.b * 0.9);
V = V + model.c * ((t - t0)^2 + (t - t0)^4);
if nargout > 1
  dF = model.c * (2*(t-t0) + 4 * (t - t0)^3);
  drho = (-model.b) * (ones(d,1) * (rho ./ s)) .* r;
  dV = dV + dF * drho;
end

end


%% TEST ROUTINE
% just test whether Vfun is implemented correctly
function test_model_toyeam()

model = model_toyeam_h(4.0, 2.0, 10, 2);
disp('-----------------------------------');
disp('    testing dV ');
disp('-----------------------------------');
test_Vfun(model);
disp('-----------------------------------');

model
model.A
model.rCB



end