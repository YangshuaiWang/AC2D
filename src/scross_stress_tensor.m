% function Sigma = scross_stress_tensor(square, U)
% 
% Input:
%	 square, T_e structure.
%        A, shear (only for testing)--> U, deformation.
% Output: Sigma, pre-computed stress tensor for T_e on elements' edges.
% Version:First release		Oct-16-2015
function [V, Sigma] = scross_stress_tensor(square, model, geom, U)

% test routine
    if nargin == 0
     test_get_stress_tensor();
     return;
    end
% check input
% if geom.geom_analyze == 0
%     error('must call geom_analyze before get_stress_tensor!');
% end
% if isfield(geom, 'volX') == 0;
%     error('mush call gqc23_prep_geom before get_stress_tensor!');
% end
% if ~strcmp(model.id, 'tri2d_toyeam_h')
%     error('consider toyeam model only by now!');
% end

% if (~isfield(geom, 'Tri'))
  geom = geom_2dtri_vacTri(geom);
% end
%% parameter set
% lattice directions, auxiliary operators
a1 = [1;0];
Q6 = [cos(pi/3), -sin(pi/3); sin(pi/3), cos(pi/3)];
a2 = Q6 * a1;
a3 = Q6 * a2;
aa = [a1, a2, a3, -a1, -a2, -a3];

dDim = 2;
nsqT = size(square.T,2);
T = zeros(dDim,dDim,nsqT);
V = zeros(1,nsqT);
% V1 = zeros(1,nsqT);

sqT = square.T;
sqX = square.X;
sqU = square.U;
X_nn = square.X_nn;
U_nn = square.U_nn;
NN = square.NN;
Tri = geom.Tri;
% g = zeros(2,6);
%% main loop: iteration over verteces.
for k = 1:size(sqX,2)
%% find interaction neighbours
    IN = find(NN(:,k));
    
%     if ismember(k, sqT(:,37371)); keyboard; end;
%     if k == 34882; keyboard; end;
    
    r = ( X_nn(:, IN) - X_nn(:, k) * ones(1,length(IN)) );
    flag = sortE(IN, r);
    u = zeros(2,6);
%     iz = find(flag == 0);
    if any(flag == 0)
%         len_iz = numel(iz);
%         fg = flag;
%         fg(fg==0) = [];
%         IN = IN(fg);
%         IN((iz(1)+len_iz):(numel(IN)+len_iz)) = IN(iz(1):numel(IN));
%         IN(iz) = 0;
        IN_tmp = zeros(6,1);
        for a = 1:6
            if flag(a) == 0;
				p_a = sqX(:,k) + aa(:,a);
                
                % BUG HERE !!!
%                 iT_a = tsearchn(geom.X', Tri, p_a');
				iT_a = geom_2dtri_aPT(p_a, geom);
                
                if isnan(iT_a)
                    u(:,a) = model.F0*aa(:,a);
                else
                    t_a = geom.T(:,iT_a);
                    x_a = geom.X(:,t_a);
                    u_a = U(:,t_a);
                    du_a = compute_stress(x_a, u_a, model);
                    u(:,a) = du_a*aa(:,a);
                end
            else
                IN_tmp(a) = IN(flag(a));
				u(:,a) = U_nn(:,IN(flag(a))) - sqU(:,k); 
            end
        end
        IN = IN_tmp;
    else
        IN = IN(flag);
        u = U_nn(:,IN) - sqU(:,k)*ones(1,length(IN));
    end
    
    %% assemble dV to elements while computing by sites.
    % 1. evaluate the pair potential part
    s = sqrt(sum(u.^2, 1));
    J = exp(-2*model.a*(s-1)) - 2 * exp(-model.a*(s-1));
%    v = 0.5*sum(J);
    dJ = -2*model.a * ( exp(-2*model.a*(s-1)) - exp(-model.a*(s-1)) );
    % forces
    dV = 0.5 * (ones(2,1) * (dJ./s)) .* u;
    % 2. evaluate the EAM part
    rho = exp(-model.b * s);
    t = sum(rho);
    t0 = model.rho0;
    v = 0.5*sum(J) + model.c * ((t - t0)^2 + (t - t0)^4);
    v = v - model.W0;

    dF = model.c * (2*(t-t0) + 4 * (t - t0)^3);
    drho = (-model.b) * (ones(2,1) * (rho ./ s)) .* u;
    dV = dV + dF * drho;
    
%   compute dV \otimes a
    sigma = zeros(dDim, dDim, length(IN));
    for i = 1:6
        sigma(:, :, i) = kron(dV(:,i),aa(:,i)');
    end

    T_k = ceil(find(sqT(:)==k)/3);
    V(T_k) = V(T_k) + v/6;    

    for i = 1:6
        j1 = flag(i);
            if j1 == 0 
                continue;
            end
        T_IN = sqT(:,T_k);
        T_IN = ceil(find(T_IN == IN(i))/3);
        Tneigs = T_k(T_IN);
        if isempty(Tneigs); continue; end;
        for j2 = 1:numel(Tneigs)
            T(:,:,Tneigs(j2)) = T(:,:,Tneigs(j2)) + sigma(:,:,i);
        end
    end
end
Sigma = 2*T/sqrt(3);
% V = V;
return;
end

%% sort edges according to lattice direction.. 
% familiar with sorted_neigs in bqc_prep_geom
% TODO: while length(r) < 6;
function flag = sortE(IN, r)
% sort the edges in site in oder (a1--a6)
Prhs = find(r(1,:)>0); Plhs = setdiff(1:length(IN), Prhs);
rrhs = r(2,Prhs); rlhs = r(2,Plhs);
if length(r) == 6
    [~, Irhs] = sort(rrhs, 2); [~, Ilhs] = sort(rlhs, 2);
    Prhs = Prhs(Irhs); Plhs = Plhs(Ilhs);
    if numel(Prhs) == 2; keyboard; end
    flag = [Prhs(2), Prhs(3), Plhs(3), Plhs(2), Plhs(1), Prhs(1)];
    return
end
flag = zeros(1,6);
for i = 1:length(rrhs)
    y = rrhs(i);
    if abs(y) < 1e-12
        flag(1) = Prhs(i);
        continue;
    end
    if y > 0
        flag(2) = Prhs(i);
    else 
        flag(6) = Prhs(i);
    end
end
for i = 1:length(rlhs)
    y = rlhs(i);
    if abs(y) < 1e-12
        flag(4) = Plhs(i);
        continue;
    end
    if y >0
        flag(3) = Plhs(i);
    else
        flag(5) = Plhs(i);
    end
end
return
end
%% test routine
function test_get_stress_tensor()
geom = geom_2dtri_mcrack(3, 5, 12, 1, 1.5);
geom = gqc23_prep_geom(geom);
figure;
geom.plot(geom);
for i = 1:geom.nT
    t = geom.T(:,i); x = sum(geom.X(1,t))/3; y = sum(geom.X(2,t))/3;
    text(x, y, num2str(i));
end
figure;
geom.plot(geom);
for n = 1:geom.nX
  text(geom.X(1,n)+0.2, geom.X(2,n), num2str(n));
end
a = 4; b = 3; c = 0; rho0 = 6 * exp(-b); rCutH = 1;
model = model_toyeam_h(a, b, c, rho0, rCutH);
A = [1 .001; .002 1];
U = A*geom.X;
Sigma = get_stress_tensor(geom, model, U);
% format long;
disp('-----------------------------------');
disp('    testing sigma ');
disp('-----------------------------------');
switch geom.fulla
    case 1
        for i = [1, 88, 107, 106]
            str = sprintf('triangulation index: %i',i);
            disp(str);
            for j = 1:6
                str = sprintf('\tdirection:  a%i ',j);
                disp(str);
                str = sprintf('\t\t%3.8f\t%3.8f\n\t\t%3.8f\t%3.8f',...
                    Sigma(1,1,j,i), Sigma(1,2,j,i), Sigma(2,1,j,i), Sigma(2,2,j,i));
                disp(str)
            end
%             disp('-----------------------------------');
            pause;
        end
    case 0
        for i = [166, 7, 122, 202, 205, 131]
            str = sprintf('triangulation index: %i',i);
            disp(str);
            str = sprintf('\t\t%3.8f\t%3.8f\n\t\t%3.8f\t%3.8f',...
                    Sigma(1,1,i), Sigma(1,2,i), Sigma(2,1,i), Sigma(2,2,i));
            disp(str);
            disp('-----------------------------------');
            pause(0.5);
        end;
    otherwise
        error('unkown value of geom.fulla');
end
end
