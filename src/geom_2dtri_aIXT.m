% function iT = geom_2dtri_aIXT(geom, iP)
% Program: find element indexed by iT with verteces' indeces within iP.
% Author : M.Liao
% Version: First released   Jan-08-2016
function iT = geom_2dtri_aIXT(geom, iX)
iT = ismember(geom.T(:),iX);
iT = find(iT == 1);
iT = ceil(iT/3);
iT = unique(iT);
iT = iT';
end