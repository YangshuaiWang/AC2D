% function Sigma = get_stress_tensor(geom, model, U)
% 
% Input: geom, geometry structure.
%        model, model structure.
%        A, shear (only for testing)--> U, deformation.
% Output: Sigma, pre-computed stress tensor.
% Version: 1.1 Keep consistent with gqc23_energy Nov-11-2014
%		   1.2 modified gqc23 part: remove reconstruction of dV, reassemble for Tic.	Oct-23-2015
%          1.3 usage of flag modified, to reduce complexity; 
%              eliminate reconstruction of dV for coupling parts; 
%              eliminate coefficient c while conbination.   Jan-13-2016
%          1.4 reconstructed continuum tensor around the interface 
%              before coupling.                  Mar-05-2016 
function Sigma = get_stress_tensor(geom, model, U)

% test routine
if nargin == 0
  test_get_stress_tensor();
  return;
end
% check input
if geom.geom_analyze == 0
    error('must call geom_analyze before get_stress_tensor!');
end
if isfield(geom, 'volX') == 0;
    error('mush call gqc23_prep_geom before get_stress_tensor!');
end
if ~strcmp(model.id, 'tri2d_toyeam_h')
    error('consider toyeam model only by now!');
end
%% parameter set
% lattice directions, auxiliary operators
a1 = [1;0];
Q6 = [cos(pi/3), -sin(pi/3); sin(pi/3), cos(pi/3)];
a2 = Q6 * a1;
a3 = Q6 * a2;
aa = [a1, a2, a3, -a1, -a2, -a3];

dDim = geom.dDim;
rDim = model.rDim;
T = zeros(dDim,dDim,geom.nT);
% Tc = zeros(dDim,dDim,geom.nT);
Tf = zeros(dDim,dDim,6,geom.nT);
%% computation of the atomistic core along with interface
% Ic = find(geom.volT > 0.5+1e-10);
% Ia = setdiff(1:geom.nT, Ic);
volX = geom.volX;
%% main loop: iteration over elements
% if geom.fulla; c_dummy(1) = c_interface(1); end; % TODO: Actually, with
% the setting of gqc23, full atom case will be considered in the same way.
Ia = find(volX~=0);
C = zeros(1,6);
for k = Ia
    %% find interaction neighbours
    % define modified modulus function
%     disp(num2str(k));
    mod6 = @(i_)(mod(i_-1, 6)+1); 
    IN = find(geom.Neigs(:,k));
    r = ( geom.X(:, IN) - geom.X(:, k) * ones(1,length(IN)) );
    flag = sortE(r);
    aET = geom.aET;
    % the column shares the same order as before. 
    [~, col] = find(geom.E == k);
    Tneigs = aET(:,col);
    iz = find(flag == 0);
    if isempty(iz)
        IN = IN(flag);
        TN = Tneigs(:,flag);
    else
        len_iz = numel(iz);
%         if len_iz > 1; keyboard; end;
        fg = flag;
        fg(fg==0) = [];
        IN = IN(fg);
        IN((iz(1)+len_iz):(numel(IN)+len_iz)) = IN(iz(1):numel(IN));
        IN(iz) = 0;
        TN = Tneigs(:, fg);
        TN(:, (iz(1)+len_iz):(size(TN,2)+len_iz)) = TN(:,iz(1):size(TN,2));
        TN(:,iz) = 0;
    end
    % cope with normal atoms and reconstruct atoms.
    switch volX(k)
        case 1
            if isempty(iz) 
                u = U(:,IN) - U(:,k)*ones(1,length(IN));
            else
                in = IN;
                in(iz) = k;
                u = U(:,in) - U(:,k)*ones(1,length(in));
            end
        case -1
            % compute the 6 C parameters
            for j = 1:6
            % if IN(j) is an interface or atomistic node then C(j) = 1.
                if abs(geom.volX(IN(j))) == 1
                    C(j) = 1;
            % otherwise C = 2/3.
                else
                    C(j) = 2/3;
                end
            end

            % compute reconstructed positions (relative to central atom)
            g = U(:, IN) - U(:, k) * ones(1,length(IN));
%             g = g(:, flag);
%             [~, id] = sort(flag);
%             C = C(flag);
            for j = 1:6
                u(:, j) = C(j) * g(:, j) ...
                        + (1-C(j)) * (g(:, mod6(j-1)) + g(:, mod6(j+1)));
            end
%             u = u(:,id);
        otherwise
            error('UNPROPER VALUE OF volX!!');
    end
    %% assemble dV to elements while computing by sites.
    % 1. evaluate the pair potential part
    s = sqrt(sum(u.^2, 1));
    dJ = -2*model.a * ( exp(-2*model.a*(s-1)) - exp(-model.a*(s-1)) );
    % forces
    dV = 0.5 * (ones(2,1) * (dJ./s)) .* u;
    itmp = isnan(dV);
    if sum(sum(itmp))~=0;
        dV(itmp) = 0;
    end
    % 2. evaluate the EAM part
    rho = exp(-model.b * s);
    t = sum(rho);
    t0 = model.rho0;

    dF = model.c * (2*(t-t0) + 4 * (t - t0)^3);
    drho = (-model.b) * (ones(2,1) * (rho ./ s)) .* u;
    itmp = isnan(drho);
    if sum(sum(itmp))~=0;
        drho(itmp) = 0;
    end
    dV = dV + dF * drho;
 
% ======================================================================= %
% This part probably is the source of all problems.
    t = zeros(size(dV));
    if volX(k) == -1;
%         dV = dV(:,flag);
        for j = 1:6
            t(:,j) = C(j) * dV(:,j) ...
                + (1-C(mod6(j-1))) * dV(:, mod6(j-1)) ...
                + (1-C(mod6(j+1))) * dV(:, mod6(j+1));
        end
%         dV = t(:,id);
        dV = t;
    end
% ======================================================================= %
%   compute dV \otimes a
    sigma = zeros(dDim, dDim, 6);
%     s_t = sigma;
    for i = 1:6
%         if IN(i) == 0 || Tneigs(1,j) == 0
%             continue;
%         end
%         if size(dV,2) < 6; 
%             disp(num2str(k));
%             keyboard; 
%         end 
        if i == iz; continue; end;
%         sigma(:, :, i) = kron(dV(:,i),aa(:,i)');
        sigma(:,:,i) = dV(:,i)*aa(:,i)';
%         s_t(:,:,i) = t(:,i)*aa(:,i)';
    end
    % TODO: check this out for small full atom core case
    %       however, those on the Bdry should be considered in another way.
    if geom.fulla
        for i = 1:6
%             j = flag(i);
%             if j == 0 || Tneigs(1,j) == 0
%                 continue;
%             end
            if i == iz; continue; end;
            if TN(1,i) ~= 0;
                Tf(:,:,i,TN(1,i)) = Tf(:,:,i,TN(1,i)) + sigma(:,:,i);
            end
            if TN(2,i) ~= 0;
                Tf(:,:,i,TN(2,i)) = Tf(:,:,i,TN(2,i)) + sigma(:,:,i);
            end
        end
%         continue;
    end

    for i = 1:6
%         j = flag(i);
%         if j == 0 || Tneigs(1,j) == 0
%             continue;
%         end
        if i == iz; continue; end;
        if TN(1,i) ~=0;
            T(:,:,TN(1,i)) = T(:,:,TN(1,i)) + sigma(:,:,i);
        end
        if TN(2,i) ~= 0;
            T(:,:,TN(2,i)) = T(:,:,TN(2,i)) + sigma(:,:,i);
        end
    end
end
if geom.fulla
    Sigma = Tf;
    return
end
%% computation of the continuum region (iteration over elements)
% loop over all elements
volT = geom.volT;
% Ta = zeros(1,0);
% for i = Ia
%     t = ceil(find(geom.T(:)==i)/3);
% %     t = t(geom.volT(t)>0.5);
%     Ta = [Ta, t'];
% end
% Ta = unique(Ta);
Ta = geom_2dtri_aIXT(geom, Ia);
% Tic = find(volT(Ta)>0 & volT(Ta)<0.5);
Tic = Ta(volT(Ta)>0 & volT(Ta)<0.5);
Tc = setdiff((1:geom.nT), Ta);
for k = Tc 
    % skip elements inside the ''interface''
    if geom.volT(k)  <= 0;
        keyboard;
    end
    % compute Du
    t = geom.T(:, k);
    J = zeros(dDim, dDim);
    Du = zeros(rDim, dDim);
    for j = 1:dDim
        J(:,j) = geom.X(:,t(j+1)) - geom.X(:,t(1));
        Du(:,j) = U(:,t(j+1)) - U(:,t(1));
    end
    if det(J) < 0.01
        error('ERROR: small element in assemble_cb');
    end
    Du = Du / J;
    % evaluate the Cauchy--Born stress tensor
    [~, dW] = model.Wfun(model, Du);
    % update gradient
    dW = reshape(dW, rDim, dDim);
    if max(T(:,:,k))~=0; keyboard; end
    T(:,:,k) = dW;
end
for k = Tic
% 	c = 2*volT(k);
    t = geom.T(:, k);
    J = zeros(dDim, dDim);
    Du = zeros(rDim, dDim);
    for j = 1:dDim
        J(:,j) = geom.X(:,t(j+1)) - geom.X(:,t(1));
        Du(:,j) = U(:,t(j+1)) - U(:,t(1));
    end
    if det(J) < 0.01
        error('ERROR: small element in assemble_cb');
    end
    Du = Du / J;
% ============================================================== %
% %  original version, seems wrong! ----M.Liao Dec 8 2015
%     [~, dW] = model.Wfun(model, Du);
%     dW = reshape(dW, rDim, dDim);
%     T(:,:,k) = T(:,:,k) + c*dW;
% ============================================================== %
	r = Du*aa;
	[~, dV] = model.Vfun(model, r);	
    
    C = 2/3*ones(1,6);
    dv = dV;
    for j = 1:6
        dv(:,j) = C(j) * dV(:,j) ...
            + (1-C(mod6(j-1))) * dV(:, mod6(j-1)) ...
            + (1-C(mod6(j+1))) * dV(:, mod6(j+1));
    end
    dV = dv;
        
	iX  = geom.volX(t);
	itmp = find(iX==0);
	r = zeros(2,0);
	for is = itmp
		xs = geom.X(:,t(is));
		ie = setdiff([1,2,3], is);
		xe = geom.X(:,t(ie));
		r = [r, xe - repmat(xs, 1,2)];
	end
	flag = sortE(r);
	idx = find(flag);
	for i = idx;
% 		T(:,:,k) = T(:,:,k) + c*dV(:,i)*aa(:,i)';
        T(:,:,k) = T(:,:,k) + dV(:,i)*aa(:,i)';
	end
end
Sigma = 2*T/sqrt(3);
return;
end

%% sort edges according to lattice direction.. 
% familiar with sorted_neigs in bqc_prep_geom
% TODO: while length(r) < 6;
function flag = sortE(r)
% sort the edges in site in oder (a1--a6)
Prhs = find(r(1,:)>0); Plhs = setdiff(1:size(r,2), Prhs);
rrhs = r(2,Prhs); rlhs = r(2,Plhs);
if length(r) == 6
    [~, Irhs] = sort(rrhs, 2); [~, Ilhs] = sort(rlhs, 2);
    Prhs = Prhs(Irhs); Plhs = Plhs(Ilhs);
    flag = [Prhs(2), Prhs(3), Plhs(3), Plhs(2), Plhs(1), Prhs(1)];
    return
end
flag = zeros(1,6);
for i = 1:length(rrhs)
    y = rrhs(i);
    if abs(y) < 1e-4
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
    if abs(y) < 1e-4
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
