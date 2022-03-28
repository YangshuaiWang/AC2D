% function geom = geom_2dtri_refinemesh(geom, IdxT)
% 
% IdxT : The indexes of triangulations in field geom.t to be refined,
%       should be a column vector.
% Author: M.Liao
% Version:
%       1.1: jigglemesh invoked. 12-Nov-2014.
%       1.2: param added. embedded function geom_refine_corners,
%       geom_refine_edges, geom_refine_jiggle written. 14-Nov-2014
%       1.21: geom_refine_corner rewrite and modify geom_refine_edge with
%       seg to construct segments correctly, then jiggle will not move the
%       outermost boundary. Mar-29-2015

function geom = geom_2dtri_refinemesh(geom, IdxT, p, v)

nargs = nargin;

if nargs == 0
  test_geom_2dtri_refinemesh();
  return;
end

if nargs == 1
    error(message('geom:refinemesh:IndexOfTrianglesNotGiven'));
end

[m, n] = size(IdxT);
% if m == 1 || n > 1;
%     error(message('geom:refinemesh:IndexOfTrianglesNotColunm'));
% end

if rem(nargs+2,2)
   error(message('geom:refinemesh:NoParamPairs'));
end

% Default value
jiggle = 'on';

if nargs > 2
    Param = eval('p');
    Value = eval('v');
    if strcmpi(Param, 'jiggle')
        jiggle = Value;
        if ~ischar(jiggle)
            error(message('geom:refinemesh:OptNotString'));
        elseif ~strcmpi(jiggle,'off') && ~strcmpi(jiggle,'on');
            error(message('geom:refinemesh:OptInvalidString'))
        end
    end  
end

%% compute corner info
% compute atoms' corner info
K0 = geom.K0;
c = zeros(1,6);
L1 = 2*K0 + 1; L2 = 1; 
c(1) = 1;
c(2) = c(1) + L2;
c(3) = c(2) + L1;
c(4) = c(3) + L2;
c(5) = c(4) + L2;
c(6) = c(5) + L1;
% c(7) = c(6) + L2 - 1;
c(7) = c(6) + L2;
% compute nodes' corner info
% Idx = geom.iBdry;
Idx = geom.iBdry;
X = geom.X(:,Idx);
uIdx = find(X(2,:)>0 | abs(X(2,:))<1e-5);
lIdx = setdiff(1:length(geom.iBdry), uIdx);
p_up = X(1,uIdx);
[~, uid] = sort(p_up, 'descend' );
p_lw = X(1,lIdx);
[~, lid] = sort(p_lw, 'ascend' );
id = [uIdx(uid), lIdx(lid)];
Idx = Idx(id);
geom.iBdry = Idx;
[C, Seg] = geom_refine_corners(geom);
%% compute decomposed geometry matrix dl;
X_A = geom.X(:,c(1:6)); X_B = geom.X(:,C(1:6));

x_A_s = X_A(1,1:6); x_A_e = [X_A(1,2:6), X_A(1,1)];
y_A_s = X_A(2,1:6); y_A_e = [X_A(2,2:6), X_A(2,1)];

x_B_s = X_B(1,1:6); x_B_e = [X_B(1,2:6), X_B(1,1)];
y_B_s = X_B(2,1:6); y_B_e = [X_B(2,2:6), X_B(2,1)];

dl_A = zeros(7,6);
dl_A(1,:) = 2;
dl_A(2,:) = x_A_s;
dl_A(3,:) = x_A_e;
dl_A(4,:) = y_A_s;
dl_A(5,:) = y_A_e;
dl_A(6,:) = 0;
dl_A(7,:) = 1;

dl_B = zeros(7,6);
dl_B(1,:) = 2;
dl_B(2,:) = x_B_s;
dl_B(3,:) = x_B_e;
dl_B(4,:) = y_B_s;
dl_B(5,:) = y_B_e;
dl_B(6,:) = 1;
dl_B(7,:) = 0;

dl = [dl_B, dl_A];
%% compute edges
e = geom_refine_edges(geom, c, C, Seg);
t = geom.T; t = [t; ones(1,size(geom.T,2))];
p = geom.X;
[p,~,t] = refinemesh(dl,p,e,t,IdxT);
%% update geom
geom.X = p;
geom.T = t(1:3, :);
PIdx = geom.nX + 1:length(geom.X);
geom.nX = length(geom.X);
geom.nT = length(geom.T);
onL = false * ones(1, geom.nX);
onL(geom.onL == true) = true;
geom.onL = onL;

PBdry = geom_Bdry_point(geom, geom.iBdry(1));
iBdry = geom_Bdry_point(geom, PIdx);
PIdx(iBdry < PBdry) = [];
% geom.iBdry = [geom.iBdry; PIdx'];
% Idx = geom.iBdry;
Idx = [geom.iBdry; PIdx'];
X = geom.X(:,Idx);
uIdx = find(X(2,:)>0 | abs(X(2,:))<1e-5);
lIdx = setdiff(1:length(Idx), uIdx);
p_up = X(1,uIdx);
[~, uid] = sort(p_up, 'descend' );
p_lw = X(1,lIdx);
[~, lid] = sort(p_lw, 'ascend' );
id = [uIdx(uid), lIdx(lid)];
geom.iBdry = Idx(id);

geom = geom_analyze(geom, 1);

if strcmpi(jiggle, 'on');
    p = geom_refine_jiggle(geom);
end
geom.X = p;
if isfield(geom, 'volT')
    geom = rmfield(geom, {'volX','volT'});
end
if isfield(geom, 'X_marker')
    X_mold = geom.X_marker;
    X_mnew = -1*ones(1,geom.nX);
    X_mnew(1:length(X_mold)) = X_mold;
    X_mnew(geom.iBdry) = 0;
    geom.X_marker = X_mnew;
end
% geom.p = p; geom.e = e; geom.t = t;
end
% =============================================================== %
%% function to compute corners of inner and outer boundaries of geom
function [CI, Seg] = geom_refine_corners(geom)
% Seg = zeros(0,1);
Idx = geom.iBdry;
% Idx = [Idx(end); Idx; Idx(1)];
% X = geom.X(:,Idx);
% % X = geom.X(:,geom.iBdry);
% uIdx = find(X(2,:)>0 | abs(X(2,:))<1e-5);
% lIdx = setdiff(1:length(geom.iBdry), uIdx);
% p_up = X(1,uIdx);
% [~, uid] = sort(p_up, 'descend' );
% p_lw = X(1,lIdx);
% [~, lid] = sort(p_lw, 'ascend' );
% id = [uIdx(uid), lIdx(lid)];
% % geom.iBdry = geom.iBdry(id);
% Idx = Idx(id);
geom.iBdry = Idx;
Idx = [Idx(end); Idx; Idx(1)];
X = geom.X(:,Idx);
detp = zeros(1,0);
% p = X(:,[1,2,3]);
% p = diff(p, [], 2);
Seg = 2;
% det_std = abs(det(p))/100;
for i = 3:size(Idx,1)-1
    p = X(:,[i-1, i, i+1]);
%     J = zeros(2, 2);
%     for j = 1:2
%         J(:, j) = p(:, j+1) - p(:, 1);
%     end
%     volT = abs(det(J)) / 2;
%     if volT < geom.N*1e-2; continue; end
    p = diff(p, [], 2);
    theta = acos(dot(p(:,1),p(:,2))/(norm(p(:,1))*norm(p(:,2))));
    if theta < pi/6; continue; end;
    detp = [detp, theta];
    Seg = [Seg; i];
end
if numel(Seg)>6; keyboard; end;
CI = Idx(Seg);
if numel(CI) < 6; keyboard; end;
n = CI(2) - CI(1);
CI(7) = CI(6) + n;
Seg = [Seg-1; length(geom.iBdry)+1];
end
% =============================================================== %
%% function to compute edges
function e = geom_refine_edges(geom, c, C, Seg)
X_A = geom.X;
X_A(:,c(7)) = X_A(:,c(1));
% X_A = X_A(:,c(1):c(7)); 
% X_A(:,c(7)) = X_A(:,c(1));
l = find(geom.onL == true,1)-1;
X_B = [geom.X, geom.X(:,C(1))];
iBdry = [geom.iBdry; geom.iBdry(1)];
e = zeros(7,0);
for i = 1:6
    seg = iBdry(Seg(i):Seg(i+1));
    x_s = X_B(:,seg(1)); x_e = X_B(:,seg(end));
    L = norm(x_e-x_s,2);
    para_e = 0;
    for j = 2:length(seg)
        x_s = X_B(:,seg(j-1)); x_e = X_B(:,seg(j));
        para_s = para_e;
        para_e = para_s + norm(x_e-x_s,2)/L;
        e = [e,[seg(j-1), seg(j), para_s, para_e, i, 1, 0]'];
    end
    if 1 - para_e >1e-5; error('LineParameterError!'); end
    e(4,end) = 1;
end
e(2,end) = C(1);

for i = 1:6
    x_s = X_A(:,c(i)); x_e = X_A(:,c(i+1));
    L = norm(x_e-x_s,2);
    para_e = 0;
    for j = c(i)+1:c(i+1)
        x_s = X_A(:,j-1); x_e = X_A(:,j);
        para_s = para_e;
        para_e = para_s + norm(x_e-x_s,2)/L;
        e = [e,[l+j-1, l+j, para_s, para_e, i+6, 0, 1]'];
    end
    if 1 - para_e >1e-5; error('LineParameterError!'); end
    e(4,end) = 1;
end
e(2,end) = c(1) + l;
end
% =============================================================== %
%% function to call jiggle
function p1 = geom_refine_jiggle(geom)
% Program: call jiggle_mesh with 'opt' and 'iter' be 'mean' and 'inf'. 
%       TODO: modify params.
%   Input: geom, proper geometry structure; p, e, t, proper geometry
%       features for pdftoolbox.
%   Output: p1, jiggled point location.
% Author: M.Liao
% Version: First release 13-Nov-2014
%       1.1. change argument to geom only. 14-Nov-2014 

% 1. dig atom core out.
m = sum(geom.onL);
c = zeros(1,6);
K0 = geom.K0; K1 = geom.K1;
c(6) = m - K1 + 1;
c(5) = c(6) - K1 - 2*K0;
c(4) = c(5) - K1;
c(3) = c(4) - K1;
c(2) = c(3) - K1 - 2*K0;
c(1) = c(2) - K1;
V = 1:(c(1) - 1);
geom_v = geom_create_vacancies(geom, V);
% 2. compute proper p, e, t as arguments of jigglemesh
[C, Seg] = geom_refine_corners(geom_v);
c = c - c(1) + 1;
c(7) = c(6) + c(2) - c(1);
e = geom_refine_edges(geom_v, c, C, Seg);
t = geom_v.T; t = [t; ones(1,size(geom_v.T,2))];
p = geom_v.X;
p1 = jigglemesh(p,e,t,'opt','mean','iter',inf);
p1 = [geom.X(:,V), p1];
end
% =============================================================== %
%% function to update iBdry.
function iBdry = geom_Bdry_point(geom, PIdx)
%% main loop
l = length(PIdx);
if l == 0;
    error('Wrong element index!');
end

X = geom.X(:,PIdx);
x = abs(X(1,:))-geom.K0; y = abs(X(2,:));
theta = angle(x + y*1i);
% 3. compute number of layers
iBdry = zeros(1,l);
for i = 1:l
    switch theta(i) < pi/3
        case 1
            iBdry(i) = floor(x(i) + tan( pi/6)*y(i)) - geom.K1;
        case 0
            iBdry(i) = floor(y(i)/sin(pi/3)) - geom.K1;
        otherwise
                error('Wrong theta Value!');
    end
end
return
end
% =============================================================== %
%% Test routine
function test_geom_2dtri_refinemesh()
geom = geom_2dtri_mcrack(3,5,20,1,1.5);
figure();
geom.plot(geom);
for i = 1:geom.nT
    t = geom.T(:,i); x = sum(geom.X(1,t))/3; y = sum(geom.X(2,t))/3;
    text(x, y, num2str(i));
end
geom = geom_2dtri_refinemesh(geom,[345 345 345 345 345]', 'jiggle', 'on');
figure();
geom.plot(geom);
end