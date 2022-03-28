% function [heart, h] = geom_2dtri_h(geom)
% Program: compute h (longest element edge) and triangle heart for each element.
% Input:
%     geom, proper geometry structure.
% Output:
%     geom.heart, triangle heart of each element.
%     geom.h, array of longest edge for every elements.
%     geom.hId, label for further usage in geom_2dtri_aPT.m
% Author: M.Liao
% Version: First release    2014 (compute h only).
%          1.10  4-Mar-2015: heart computed. 
function geom = geom_2dtri_h(geom)
h = zeros(1,geom.nT);
e = zeros(3,geom.nT);
X = geom.X;
T = geom.T;
for k = 1:geom.nT
    t = T(:,k);
    x = X(:,t);
    e_k = [x(:,3) - x(:,2), x(:,3) - x(:,1), x(:,2) - x(:,1)];
    e_k = sqrt(sum(e_k.^2))';
    h(k) = max(e_k);
    e(:,k) = e_k;
end
heart = geom_triangle_heart(X, T);
geom.heart = heart;
geom.maxH = h;
geom.hId = size(heart,2);
geom.LenE = e;
return;
end

%% function x_h = geom_triangle_heart(k, geom)
% compute coordinate of triangle heart x_h of triangle with index k on geom.
function heart = geom_triangle_heart(X, T)
heart = zeros(2,size(T,2));
for k = 1:size(T,2)
    x_t = X(:,T(:,k));
    r = diff([x_t, x_t(:,1)],1,2);
    r = sqrt(sum(r.^2,1)); % r = [c a b]
    r = [r([2,3,1]); r([2,3,1])];
    x_h = r.*x_t/sum(r(1,:));
    heart(:,k) = sum(x_h,2);
end
end