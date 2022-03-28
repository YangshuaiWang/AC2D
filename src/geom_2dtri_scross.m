% function [I_cross, OverlapArea] = geom_2dtri_scross(geom)
% Program: Find triangles on full atom config locating on elements'
%       boundaries. triangles save in pairs as squares under cartesian coordinates.
%     Input: geom, proper geometry structure.
%     Output: I_cross, the top-right coordinates of square contains to
%       triangles; OverlapArea, overlapping areas of T_e and T_h.
% Author: M.Liao
% Version: First release    Nov-25-2014    (isolated from error_model.m)
%     1.01 replace output with structure squre, contains idx, and oa.
%       Nov-26-2014.
function square = geom_2dtri_scross(geom)
% TODO: test routine.
if nargin == 0;
    test_geom_2dtri_scorss();
    return;
end
% TODO: check inputs (geom.label?)
if (~isfield(geom, 'volT'))
  error('ERROR: must call bqc_prep_geom before get_model_error!');
end

OverlapArea = zeros(3,0);
I_cross = zeros(2,0);
iAtom = find(geom.volX ~= 0);
iTa = geom_2dtri_aIXT(geom, iAtom);
iTc = setdiff(1:geom.nT, iTa);
% Ic = find(geom.volT > 0.5+1e-10);
detA = det(geom.A);
for k = iTc
%     if k == 1917;
%         disp(num2str(k));
%     end
    t = geom.T(:,k);
    X = geom.X(:,t);
    X = geom.A\X;
    x = X(1,:); y = X(2,:);

    % find effective squares
    xt = floor(min(x)); yt = floor(min(y));
    x = x-xt; y = y-yt;

    I3 = get_square([x([1,2]);y([1,2])]);
    I2 = get_square([x([1,3]);y([1,3])]);
    I1 = get_square([x([3,2]);y([3,2])]);
    I = unique([I1, I2, I3]', 'rows')';
    % compute overlapping area
    x = x(end:-1:1); y = y(end:-1:1);
    for idx = 1:size(I,2)
        i = I(1,idx); j = I(2,idx);
        x_up = [i-1, i, i]; y_up = [j, j, j-1];
        a_up = get_overlapping_area(x, y, x_up, y_up);
        x_dn = [i, i-1, i-1]; y_dn = [j-1, j-1, j];
        a_dn = get_overlapping_area(x, y, x_dn, y_dn);
        OverlapArea = [OverlapArea, [k; detA*a_up; detA*a_dn]];
    end
    I = I + [xt*ones(1,size(I,2)); yt*ones(1,size(I,2))];
    I_cross = [I_cross, I]; 
end
[I_cross, ~, map_I2OA] = unique(I_cross', 'rows');
I_cross = I_cross';
% fields of square. 
% 1. top-right point index.
square.I = I_cross;
% 2. overlapping area to T_h.
square.oa = OverlapArea;
% 3. map of square to corresponding overlapping area
square.aIoa = map_I2OA;
% 4. id related to corresponding geom.
square.id = geom.nX;
end
%% test routine
function test_geom_2dtri_scorss()
geom = geom_2dtri_mcrack(3, 3, 40, 1, 1.5);
geom = geom_analyze(geom, 1);
geom = gqc23_prep_geom(geom);
geom.plot(geom);
hold on;
tic;
square = geom_2dtri_scross(geom);
toc;
I = square.I;
for idx = 1:length(I)
    i = I(1,idx); j = I(2,idx);
    xf = [i-1,i,i,i-1]; yf = [j-1,j-1,j,j];
    X = geom.A*[xf; yf];
    fill(X(1,:), X(2,:), 'b');
end
hold off;
end
