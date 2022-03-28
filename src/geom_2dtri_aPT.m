% function aPT = geom_2dtri_aPT(X, geom)
% Program: for input coordinates X, find corresponding triangle in geom
%       that contains X.
%   Output: elements' indexes or NaN for that not inside the geometry. 
% Author: M.Liao
% Version: First release    Nov-30-2014
%   1.01, 4-Mar-2015: a box around X added to save computation berden. 
function aPT = geom_2dtri_aPT(X, geom)
% TODO: test routine.
% check input
if isfield(geom, 'hId')
   if geom.hId ~= geom.nT
       warning('geom:hId','hId changed!');
       geom = geom_2dtri_h(geom);
   end
else
    geom = geom_2dtri_h(geom);
end
% TODO: size of X should be (2, n).
n = size(X,2);
aPT = zeros(1,n);
T = geom.T;
X_T = geom.X;
heart = geom.heart;
maxH = max(geom.maxH);

for i = 1:n
    p = X(:,i);
%     p_box = [left, right, bottom, top] margins;
    p_box = [p(1)-maxH, p(1)+maxH, p(2)-maxH, p(2)+maxH];
    p_tmp = find(p_box(1)<heart(1,:) & heart(1,:)<p_box(2));
    p_tmp = p_tmp(p_box(3)<heart(2,p_tmp) & heart(2,p_tmp)<p_box(4));
    
%     r = repmat(X(:,i),1,length(p_tmp)) - heart(:,p_tmp);
%     r = bsxfun(@minus, X(:,i), heart(:,p_tmp));
%     r = sum(r.^2);
    r = pdist2(X(:,i)', heart(:,p_tmp)');
    [~, I] = sort(r);
    for j = 1:length(p_tmp)
        t = T(:,p_tmp(I(j)));
        X_t = X_T(:,t);
        if isinsidetriangle(X(:,i), X_t)
            aPT(i) = p_tmp(I(j));
            break
        end
    end
%     t = tsearchn(geom.X', geom.T', p');
%     if t ~= aPT(i)&&aPT(i) ~= 0; keyboard; end
    if aPT(i) == 0; aPT(i) = NaN; end
end

return
end
%% function output = isinsidetriangle(x, t, X)
% check if point x inside triangle whose vertexes stored in X or not.
function output = isinsidetriangle(x, X)
    V0 = (X(:,3) - X(:,1))';
    V1 = (X(:,2) - X(:,1))';
    V2 = (x - X(:,1))';
    Demo = (V0*V0')*(V1*V1') - (V0*V1')*(V1*V0');
    u = ((V1*V1')*(V2*V0') - (V1*V0')*(V2*V1'))/Demo;
    if (abs(u) < 1e-10); u = 0; end;
    if (abs(u-1) < 1e-10); u = 1; end;
    if (u < 0 || u > 1)
        output = false;
        return
    end
    v = ((V0*V0')*(V2*V1') - (V0*V1')*(V2*V0'))/Demo;
    if (abs(v) < 1e-10); v = 0; end
    if (abs(v-1) < 1e-10); v = 1; end;
    if (v < 0 || v > 1)
        output = false;
        return
    end    
    output = (u + v <= 1+1e-10);
end