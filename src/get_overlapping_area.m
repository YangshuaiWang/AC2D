% function a = compute_ovelapping_area(x1, y1, x2, y2)
% computing overlapping area of two triangles by applying polybool and polyarea
%   input: x1, y1 saves coordinates of one triangle, x2, y2 saves the others.
%   output: overlapping area
function a = get_overlapping_area(x1, y1, x2, y2)
%     this part is caused by the order of points in geom.T
%     x1 = x1(end:-1:1); y1 = y1(end:-1:1);
    [xi, yi] = polybool('intersection', x1, y1, x2, y2);
    a = polyarea(xi, yi);
    return;
end