% function geom = geom_2dtri_vacTri(geom)
% program: construct Tri as a field of geom to apply tsearchn properly.
function geom = geom_2dtri_vacTri(geom)

% in geom_analyze.m file (for edge)
ix = ceil(find(geom.aET==0)/2);
ix = geom.E(:,ix);
ix = unique(ix);
ix = setdiff(ix, geom.iBdry);

x = geom.X(:, ix)';
tri = delaunayn(x);

Tri = [geom.T'; ix(tri)];
geom.Tri = Tri;
end