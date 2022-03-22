% function iN = sorted_neigs_2dtri(i0, geom)
%
% If the geometry is the 2D triangular lattice, then this routine
% computes the six neighbours of the vertex i0, in counterclockwise
% order. If there are not 6 neighours, then in all likelihood there is
% a bug in the routine from which sorted_neigs_2dtri is called, and
% hence an error message is generated.
%
% Input: 
%      i0 : index of central node in geom.X
%    geom : valid geometry structure, associated with the 2dtri lattice
% Output:
%      iN : list of 6 indices in geom.X, which are the 6 neighbours of
%           sorted in counterclockwise order
%


function iN = sorted_neigs_2dtri(i0, geom)
iN = find(geom.NN(:, i0));
if length(iN) ~= 6
  error('ERROR: sorted_neigs_2dtri assumes there are 6 neighbours!');
end

x = geom.X(:, iN) - geom.X(:,i0) * ones(1,6);
x1 = x(:,1); x1p = [-x1(2); x1(1)];

ii = zeros(1,6); ii(1) = 1;
for i = 2:6
  if x1'*x(:,i) > 0
    if x1p'*x(:,i) > 0
      ii(2) = i;
    else
      ii(6) = i;
    end
  else
    b = x1p'*x(:,i);
    if abs(b) < 1e-5
      ii(4) = i;
    elseif b > 0
      ii(3) = i;
    else 
      ii(5) = i;
    end
  end
end
iN = iN(ii);
end
