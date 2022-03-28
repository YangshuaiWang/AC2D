% function square = geom_interp_scorss(geom, square, U, model)
% Program: compute deformed coordinates of T_e crossing elements' boundary.
%     Input: geom, proper geometry structure; I_cross, coordinates of
%       top-right vertex of squares containing T_e; U, deformed nodal
%       coordinates.
%     Output: U_Te, deformed coordinates of T_e.
% Author: M.Liao
% Version: First release    Nov-25-2014    (isolated from error_model.m)
%     1.1 replace output by structure square, fields U, X, T, aTI added
%       Nov-27-2014.
%     1.2 modify computation of U_pt. Oct-15-2015. 
function square = geom_interp_scross(geom, square, U, model)
% TODO: test routine.
if nargin == 0;
    test_geom_interp_scross();
    return;
end
% TODO: check inputs (geom.label?)
if (~isfield(geom, 'volT'))
  error('ERROR: must call bqc_prep_geom before get_model_error!');
end

% if (~isfield(geom, 'Tri'))
  geom = geom_2dtri_vacTri(geom);
% end

if isfield(geom, 'hId')
   if geom.hId ~= geom.nT
       warning('geom:hId','hId changed!');
       geom = geom_2dtri_h(geom);
   end
else
    geom = geom_2dtri_h(geom);
end 
Du = zeros(4,geom.nT);
for i = 1:geom.nT
    t = geom.T(:,i);
    x = geom.X(:,t);
    u = U(:,t);
    du = compute_stress(x, u, model);
    Du(:,i) = reshape(du, 4, 1);
end

I_cross = square.I;
OA = square.oa;
aIOA = square.aIoa;

I_min = min(min(I_cross)) - 2; I_max = max(max(I_cross));
map_s2t = spalloc(I_max-I_min, I_max-I_min, 10*(I_max-I_min));
U_Te = zeros(2,0);
X_Te = zeros(2,0);
T_e = zeros(3,0);
aEI = zeros(2,0);
Ind = zeros(1,4);
Tri = geom.Tri;
for i = 1:length(I_cross)
    % skip square include by numerical reason.
%     if i == 365; keyboard; end;
    oa = OA(2:3,aIOA == i);
    if sum(sum(oa)) < 1e-6;
        continue;
    end
    
    r = I_cross(1, i); c = I_cross(2, i);
    I_sqr = [[r-1, r, r, r-1]; [c, c, c-1, c-1]];
    X_sqr = geom.A*I_sqr;
    I_sqr = I_sqr - I_min;
    
%     T_lcted = geom_2dtri_aPT(X_sqr, geom);
    T_lcted = tsearchn(geom.X', Tri, X_sqr');
    if sum(isnan(T_lcted))~=0
%         T_tmp = geom_2dtri_aPT(X_sqr, geom);
%         if sum(isnan(T_tmp))~=0
%             continue;
%         else
%             T_lcted = T_tmp;
%         end
        continue;
    end

    tmp = zeros(1,0);
    for j = 1:4
        idx = full(map_s2t(I_sqr(1,j), I_sqr(2,j)));
        if idx ~= 0; 
            Ind(j) = idx;
            tmp = [tmp, j];
            continue; 
        end;
        X_Th = geom.X(:,geom.T(:,T_lcted(j)));
        U_Th = U(:,geom.T(:,T_lcted(j)));
        du = reshape(Du(:,T_lcted(j)),2,2);
        X_pt = X_sqr(:,j);
        
%         u_th = U(:,geom.T(:,T_lcted(j)));
%         if norm(U_Th - u_th,inf) > 1e-10;
%             keyboard;
%         end
%         J  = [X_Th(:,2)-X_Th(:,1), X_Th(:,3)-X_Th(:,1)];
%         J1 = [X_Th(:,2)-X_pt, X_Th(:,3)-X_pt];
%         J2 = [X_Th(:,3)-X_pt, X_Th(:,1)-X_pt];
%         J3 = [X_Th(:,1)-X_pt, X_Th(:,2)-X_pt];
%         U_pt = (abs(det(J1))*U_Th(:,1) + abs(det(J2))*U_Th(:,2) +...
%             abs(det(J3))*U_Th(:,3))/abs(det(J));
        U_pt = U_Th(:,1) + du*(X_pt-X_Th(:,1));
        U_Te = [U_Te, U_pt];
        map_s2t(I_sqr(1,j), I_sqr(2,j)) = size((U_Te),2);
        Ind(j) = size(U_Te,2);
    end

    if ~isempty(tmp); X_sqr(:,tmp) = []; end;
    X_Te = [X_Te, X_sqr];
%     T_e = [T_e, Ind(1:3)', Ind([1,3,4])'];
    T_e = [T_e, Ind([3,1,4])', Ind([3,2,1])'];
% %     =============================================================      %
% % for debug
%     x = U_Te(:,Ind([3,1,4]));
%     e = [x(:,3) - x(:,2), x(:,3) - x(:,1), x(:,2) - x(:,1)];
%     h = max(sqrt(sum(e.^2)));
%     if h > 2; keyboard; end
% %     =============================================================      %
    aEI = [aEI, [i;3], [i;2]];
end
aXT = tsearchn(geom.X', Tri, X_Te');
% aXT = geom_2dtri_aPT(X_Te, geom); 
% update square and return
square.U = U_Te;
square.X = X_Te;
square.T = T_e;
square.aTI = aEI;
square.aXT = aXT';
end
%% test routine
function test_geom_interp_scross()
geom = geom_2dtri_mcrack(3, 3, 20, 1, 1.5);
geom = geom_analyze(geom, 1);
geom = gqc23_prep_geom(geom);
geom = geom_2dtri_h(geom);
model = model_toyeam_h(4, 3, 0, 6*exp(-3), 1);
A = [1 .001; .002 1];
U = A*geom.X;
sqr = geom_2dtri_scross(geom);
tic;
sqr = geom_interp_scross(geom, sqr, U, model);
toc;
pub_vis_alex(geom, U, 20, 0, 0, 1.2)
hold on
triplot(sqr.T', sqr.U(1,:)', sqr.U(2,:)', 'g');
hold off
end
