% function square = geom_scross_neigs(geom, square, U, model)
% Program: After calling geom_2dtri_scross & geom_interp_scross, this code 
%     construct new fields X_nn saves the info of nearest neighbour and 
%     U_nn stores displacements of cooresponding X_nn.
%     Caution: geom.K1 >= 4
% Author : M.Liao
% Version: First release 	Oct-15-2015
% 		1.01 Enlarge Teic to include Tia, save as TeI.	Oct-21-2015
%		1.02 new computation of NN. 					Oct-22-2015
%       1.2 modified to improve performance             Mar-16-2016
function square = geom_scross_neigs(geom, square, U, model)
% check input
if ~(isfield(square, 'U'))
	warning('MATLAB:AUGMENTS', 'call geom_interp_scross first!!');
	square = geom_interp_scross(geom, square, U, model);
end

% if (~isfield(geom, 'Tri'))
  geom = geom_2dtri_vacTri(geom);
% end

% check input
if isfield(geom, 'hId')
   if geom.hId ~= geom.nT
       warning('geom:hId','hId changed!');
       geom = geom_2dtri_h(geom);
   end
else
    geom = geom_2dtri_h(geom);
end

% sqU = square.U;
X = geom.X;
T = geom.T;
sqX = square.X;
sqU = square.U;
sqT = square.T';
aXT = square.aXT;
Tri = geom.Tri;

NN = sparse( [sqT(:,1); sqT(:,1); sqT(:,2); sqT(:,2); sqT(:,3); sqT(:,3)], ...
             [sqT(:,2); sqT(:,3); sqT(:,1); sqT(:,3); sqT(:,1); sqT(:,2)], ...
             0.1*ones(6*size(sqT,1), 1) );
NN = ceil(NN);

geomVolX = find(geom.volX==-1);
X_intf_I = X(:,geomVolX);
lGeomVolX = numel(geomVolX);
Idx = zeros(1,lGeomVolX);
sqx = zeros(2,lGeomVolX);
squ = zeros(2,lGeomVolX);
axt = zeros(1,lGeomVolX);
lsqx = size(sqX, 2);
n = 0;
for i = 1:length(X_intf_I)
	tn = ceil(find(T(:)==geomVolX(i))/3);
    ixtmp = intersect(tn', aXT);
    if size(ixtmp,1)>1; keyboard; end
    isqX = zeros(1,0);
    
    for j = ixtmp
        if numel(j)~=1; keyboard; end;
        isqX = [isqX, find(aXT == j)];
    end
    isqX = unique(isqX);
    tmp = pdist2(X_intf_I(:,i)', sqX(:,isqX)');%,'euclidean', 'Smallest', 1);
    idx = find(tmp<1e-10);
    if numel(idx)>1;
        keyboard;
    end
    if isempty(idx);
        n = n + 1;
        sqx(:,n) = X_intf_I(:, i);
        squ(:,n) = U(:,geomVolX(i));
        idx = n + lsqx;
        axt(n) = -1;
    else
        idx = isqX(idx);
        aXT(idx) = -1;
    end
	
    Idx(i) = idx;
end
if n ~= 0;
    tmp = sum(sqx);
    sqx(:,tmp==0) = [];
    sqX = [sqX, sqx];
    squ(:,tmp==0) = [];
    sqU = [sqU, squ];
    axt(tmp==0) = [];
    aXT = [aXT, axt];
end
map = [geomVolX; Idx];

geomdi0 = find(geom.di==0);
X_di_0 = X(:,geomdi0);
Idx = zeros(1,length(X_di_0));
for i = 1:size(X_di_0,2)
	tn = ceil(find(T(:)==geomdi0(i))/3);
	ixtmp = intersect(tn', aXT);
    isqX = zeros(1,0);
    for j = ixtmp
        isqX = [isqX, find(aXT==j)];
    end
    isqX = unique(isqX);
    tmp = pdist2(X_di_0(:,i)', sqX(:,isqX)');%,'euclidean', 'Smallest', 1);
    idx = find(tmp<1e-10);
    if numel(idx) ~= 1;
		keyboard;
    end
 	idx = isqX(idx);
	aXT(idx) = 0;
    Idx(i) = idx;
end
map = [map,[geomdi0; Idx]];

iXia = find(geom.di==2);
ltmp = size(sqX,2);
sqX = [sqX, X(:,iXia)];
sqU = [sqU, U(:,iXia)];
aXT = [aXT,(-2)*ones(1, numel(iXia))];

ltmp = (ltmp+1):size(sqX,2);
map = [map, [iXia;ltmp]];

Du = zeros(4,geom.nT);
for i = 1:geom.nT
    t = geom.T(:,i);
    x = geom.X(:,t);
    u = U(:,t);
    du = compute_stress(x, u, model);
    Du(:,i) = reshape(du, 4, 1);
end

X_nn = [sqX, X(:,geom.di==3)];
U_nn = [sqU, U(:,geom.di==3)];

lNN = size(X_nn,2);
row = zeros(1,3*size(X_nn,2));
col = zeros(1,3*size(X_nn,2));
% computation domain decomposition:
iXeC = find(aXT>0);
iXe0 = find(aXT==0);
iXe1 = find(aXT==-1);
iXe2 = find(aXT==-2);
iXe3 = (size(sqX,2)+1):size(X_nn,2);
x_nn = zeros(2,6*numel(iXeC));
u_nn = zeros(2,6*numel(iXeC));
axt = zeros(2,6*numel(iXeC));
i_nn = 0;
i_n = 0;
for i = iXeC;
	x = sqX(:,i);
    nn = repmat(x,1,6) + geom.aa(:,1:6);
	in = find(NN(:,i));
    if numel(in) == 6; continue; end;
	r = X_nn(:,in) - repmat(x,1,numel(in));
	flag = sortE(r);
	ordn = find(flag==0);
    for j = ordn
        iTh = tsearchn(geom.X', Tri, nn(:,j)');
%         iTh = iTh';
% 		iTh = geom_2dtri_aPT(nn(:,j), geom);
        if isnan(iTh); continue; end;
        rTh = find(aXT==iTh);
        tmp = pdist2([X_nn(:,rTh),x_nn(:,1:i_nn)]', nn(:,j)');%,'euclidean', 'Smallest', 1);
        idx = find(tmp<1e-10);
        if isempty(idx);
            i_nn = i_nn + 1;
            if i_nn > size(x_nn,2); keyboard; end;
            x_nn(:,i_nn) = nn(:,j);
			X_Th = X(:,T(:,iTh));
			U_Th = U(:,T(:,iTh));
       		du = reshape(Du(:,iTh),2,2);
			U_pt = U_Th(:,1) + du*(nn(:,j)-X_Th(:,1));
            u_nn(:,i_nn) = U_pt;
			axt(i_nn) = iTh;
			idx = lNN + i_nn;
        else
            if numel(rTh) < idx;
                idx = lNN + idx - numel(rTh);
            else
                idx = rTh(idx);
            end
        end
        if i_n > size(row,2); keyboard; end;
        i_n = i_n + 1;
        row(i_n) = idx;
		col(i_n) = i;		
    end
end
if i_nn ~= 0;
    tmp = sum(x_nn);
    x_nn(:,tmp==0) = [];
    X_nn = [X_nn, x_nn];
    u_nn(:,tmp==0) = [];
    U_nn = [U_nn, u_nn];
end
iXe_01 = [iXe0, iXe1];
X_tmp = X_nn(:,iXe_01);

lNN = size(X_nn,2);

x_nn = zeros(2,6*numel(iXe0));
u_nn = zeros(2,6*numel(iXe0));
i_nn = 0;
for i = iXe0
	x = sqX(:,i);
    nn = repmat(x,1,6) + geom.aa(:,1:6);
	
	in = find(NN(:,i));
    if numel(in) == 6; continue; end;
	r = X_nn(:,in) - repmat(x,1,numel(in));
	flag = sortE(r);
	ordn = find(flag==0);
    for j = ordn
        tmp = pdist2(X_tmp', nn(:,j)');%,'euclidean', 'Smallest', 1);
        idx = find(tmp<1e-10);
        if isempty(idx);
			keyboard;
            T_lcted = tsearchn(geom.X', Tri, nn(:,j)');
% 			T_lcted = geom_2dtri_aPT(nn(:,j), geom);
            
            if isnan(T_lcted); continue; end;
            i_nn = i_nn + 1;
            x_nn(:,i_nn) = nn(:,j);
			X_Th = X(:,T(:,T_lcted));
			U_Th = U(:,T(:,T_lcted));
       		du = reshape(Du(:,T_lcted),2,2);
			U_pt = U_Th(:,1) + du*(nn(:,j)-X_Th(:,1));
            u_nn(:,i_nn) = U_pt;
			idx = lNN + i_nn;
        else
            idx = iXe_01(idx);
        end
        i_n = i_n + 1;
        row(i_n) = idx;
        col(i_n) = i;
    end
end
if i_nn ~= 0;
    tmp = sum(x_nn);
    x_nn(:,tmp==0) = [];
    X_nn = [X_nn, x_nn];
    u_nn(:,tmp==0) = [];
    U_nn = [U_nn, u_nn];
end
iXe_012 = [iXe_01, iXe2];
X_tmp = X_nn(:,iXe_012);
sNN = size(NN,2);
for i = iXe1
	x = sqX(:,i);
    nn = repmat(x,1,6) + geom.aa(:,1:6);
	
    if i > sNN; 
        ordn = 1:6;
    else
        in = find(NN(:,i));
        if numel(in) == 6; continue; end;
        r = X_nn(:,in) - repmat(x,1,numel(in));
        flag = sortE(r);
        ordn = find(flag==0);
    end
    for j = ordn
        tmp = pdist2(X_tmp', nn(:,j)');%,'euclidean', 'Smallest', 1);
        idx = find(tmp<1e-10);
        if isempty(idx);
			keyboard;
        end
        idx = iXe_012(idx);
        i_n = i_n + 1;
        row(i_n) = idx;
        col(i_n) = i;
    end
end

%%% YS: BUG here!
iXe_123 = [iXe1, iXe2, iXe3];
X_tmp = X_nn(:,iXe_123);
x_nn = zeros(2,6*numel(iXe2));
u_nn = zeros(2,6*numel(iXe2));
i_nn = 0;
lNN = size(X_nn, 2);
for i = iXe2
	x = sqX(:,i);
    nn = repmat(x,1,6) + geom.aa(:,1:6);
	
    for j = 1:6
		tmp = bsxfun(@minus, X_tmp, nn(:,j));
		idx = find(all(abs(tmp)<1e-5), 1);
        if isempty(idx);
%             BUG here
%             keyboard;
%             T_lcted = tsearchn(geom.X', Tri, nn(:,j)');
			T_lcted = geom_2dtri_aPT(nn(:,j), geom);
            
            if isnan(T_lcted); continue; end;
            
            i_nn = i_nn + 1;
            
             x_nn(:,i_nn) = nn(:,j);
			X_Th = X(:,T(:,T_lcted));
 			U_Th = U(:,T(:,T_lcted));
       		du = reshape(Du(:,T_lcted),2,2);
			U_pt = U_Th(:,1) + du*(nn(:,j)-X_Th(:,1));
            u_nn(:,i_nn) = U_pt;
            
			
            idx = lNN + i_nn;
        else
            idx = iXe_123(idx);
        end
        i_n = i_n + 1;
        row(i_n) = idx;
        col(i_n) = i;
    end
end
row(row==0) = [];
col(col==0) = [];

[ir, ic] = ind2sub(size(NN), find(NN));
row = [ir', row];
col = [ic', col];
NN = sparse(row, col, 1);

itmp = geom_2dtri_aIXT(geom, geomVolX);

tmp = geom.T(:,itmp);
Map = sparse(map(1,:), map(2,:), 1);
for i = 1:numel(tmp)
    X_idx = find(Map(tmp(i),:));
    tmp(i) = X_idx;
end

tmp = sort(tmp);
T = sort(square.T);
t = setdiff(tmp', T', 'rows');
square.X = sqX;
square.U = sqU;
square.X_nn = X_nn;
square.U_nn = U_nn;
square.NN = NN;
square.T = [square.T, t'];
[~,iTI,iTeI] = intersect(tmp', [T,t']', 'rows');
square.iTI = itmp(iTI)';
square.iTeI = iTeI';
return
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