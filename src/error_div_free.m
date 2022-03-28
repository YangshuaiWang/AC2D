% function [ufs, divFC] = error_div_free(geom, square, STh, STe)
% Program: For input geometrical structures, and STh (stress of continuum)
%       and STe (stress of atomistic), compute ufs (unit face-vector structure)
%       and divFC (divergence free constants).
% Author : M.Liao
% Version: First release    Dec-9-2015
function [S, ufs, divFC] = error_div_free(geom, square, STh, STe)
% TODO: check inputs and test_routine
% test routine
if nargin == 0
  test_error_div_free();
  return;
end
% 1. find corresponding elements, and compute stress difference.
% iTI = square.iTI;
% iTIfc = iTI(geom.volT(iTI)>0);
iNI = find(geom.volX == -1);
iTIfc = geom_2dtri_aIXT(geom, iNI);

lTIfc = numel(iTIfc);
S = zeros(2,2,lTIfc);
for i = 1:lTIfc
    it = iTIfc(i);
%     id = find(square.iTI == it);
    ite = square.iTeI(square.iTI == it);
    S(:,:,i) = STe(:,:,ite) - STh(:,:,it);
end

% 2. find related edges, sort them, find effective edges; count them based 
% on number of effective elements they lying, sort again for further use.
iFIfc = geom.aTE(:,iTIfc);
iFtmp = unique(iFIfc(:));

edge = geom.E;
volX = geom.volX;
idx = zeros(1,0);

for i = 1:numel(iFtmp)
    ie = iFtmp(i);
    e = edge(:,ie);
    if all(volX(e) == 0);
        idx = [idx, find(iFIfc == ie)];
    end
end

iFIfc(idx) = 0;
[iw, ~, in] = unique(iFIfc(:));

% 3. main loop: assembling stiffness matrix and r.h.s vector.
ROW = zeros(1,0);
COL = zeros(1,0);
row = 0;
A = zeros(1,0);
b = zeros(0,1);

ufs = zeros(4,0);

for i = 1:size(iFIfc,2)
    it = iTIfc(i);
    s = S(:,:,i);
    iT = geom.T(:,it);
    idx = (i-1)*3 + [1,2,3];
    f_idx = in(idx)';
    inull = find(f_idx == 1);

    switch numel(inull)
        case 0
            t = zeros(3,2);
            for l = 1:3
                ix = setxor([1,2,3], l);
                if isequal(ix, [1,3])
                    ix = [3,1];
                end
                x = geom.X(:,iT(ix));
                x = diff(x,[],2);
                x = x/sqrt(x'*x);
                t(l,:) = x;
                ufs = [ufs, [it; iw(f_idx(l)); x]];
            end
        case 1
            t = zeros(2,2);
            m = 0;
            for l = 1:3
                if l == inull; continue; end;
                m = m + 1;
                ix = setxor([1,2,3], l);
                if isequal(ix, [1,3])
                    ix = [3,1];
                end
                x = geom.X(:,iT(ix));
                x = diff(x,[],2);
                x = x/sqrt(x'*x);
                t(m,:) = x;
                ufs = [ufs, [it; iw(f_idx(l)); x]];
            end
        otherwise
            keyboard;
    end
    f_idx(inull) = [];
    for p = 1:numel(f_idx)
        for q = 1:2
            row = row + 1;
%             ROW = [ROW, row];
            for k = 1:numel(f_idx)
                ROW = [ROW, row];
                col = 2*(f_idx(k)-2) + q;
                COL = [COL, col];
                A_val = 2*t(p,:)*t(k,:)';
                A = [A, A_val];
            end
            b_val = 2*s(q,:)*t(p,:)';
            b = [b; b_val];
%             b(row) = b_val;
        end
    end
end


% 4. Solve the least square AC = b; and construct output.
SA = sparse(ROW, COL, A);
C = SA\b;
% C = A\b;
C = reshape(C, 2, []);
divFC = [iw(2:end)'; C];
S = reshape(S, 4, lTIfc);
S = [iTIfc; S];
return
end
%% test routine
function test_error_div_free()
load Mar7
% ex_mcrack;
% STh = get_stress_tensor(geom, model, Y);
[S, ufs, divFC] = error_div_uniform(geom, STh, STh(:,:,1));
eta = zeros(4,size(S,2));
for i = 1:size(S,2)
    s = S(2:5,i);
    iT = S(1,i);
    it = find(ufs(1,:) == iT);
    ie = ufs(2,it);
    for j = 1:numel(ie)
        ic(j) = find(divFC(1,:)==ie(j));
    end
    dfc = zeros(2,2);
    for j = 1:numel(ie)
        dfc = dfc + divFC(2:3,ic(j))*ufs(3:4, it(j))';
    end
    eta(:,i) = reshape(dfc, 4, 1) - s;
end
fprintf('----------------------------------------------------------------------------\n');
fprintf('The maximum of difference of norm is %3.5g (Uniform deformation).\n', max(sqrt(sum(eta.^2))));
end