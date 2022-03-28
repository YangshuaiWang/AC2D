
% geom = geom_create_vacancies(geom, V)
%
% removes the vertices with indices V from the mesh
%
% Input
%   geom : valid geometry structure
%      V : integer array, 1 <= V(i) <= geom.nX, of indiced to be removed
%          from geom.X
% Output
%   geom : geometry structure with vacancies removed. If geom_analyze was 
%          run on geom before, then it will automatically be called again

function geom = geom_create_vacancies(geom, V)

if nargin == 0
  test_geom_create_vacancies();
  return
end

% check input
nX = size(geom.X, 2); nT = size(geom.T, 2);
V = sort(unique(V));
if (V ~= round(V))
  error('ERROR: geom_create_vacancies requires V to be an integer array');
end

if ( (V(1) < 1) || (V(end) > nX) )
  error('ERROR: geom_create_vacancies requires 1 <= V(i) <= nX');
end



%% 1. fix the node array
% create arrays to store the new nodes and node indices
nV = length(V);
I = 1:nX;
X = zeros(2, nX - nV);
% add an extra index to make the subsequent loop work
V = [V(:)', nX+1];
% loop through vacancies
iX = V(1)-1;
X(:, 1:iX) = geom.X(:, 1:iX);
for jV = 1:nV
  I(V(jV)) = 0;
  oldinds = (V(jV)+1):(V(jV+1)-1);
  newinds = (iX+1):(iX+V(jV+1)-V(jV)-1);
  I(oldinds) = newinds;
  
  
   if isnan(geom.X(:, oldinds))
        keyboard;
   end
    
  X(:, newinds) = geom.X(:, oldinds);
  iX = iX + length(oldinds);
end

%% 2. fix the element array
T = zeros(size(geom.T));
iTnew = 0;
for iTold = 1:nT
  % if element needs to be removed do nothing. if it is kept, add it to T
  if min( I( geom.T(:, iTold) ) ) > 0
    iTnew = iTnew + 1;
    T(:, iTnew) = I( geom.T(:, iTold) );
  end
end

%% 3. fix geom structure
geom.X = X; geom.nX = size(geom.X, 2);
geom.T = T(:, 1:iTnew);  geom.nT = size(geom.T, 2);
geom.onL(I == 0) = [];
% fix iPer and iBdry
% WARNING: this assumes that no boundary vertex was removed!
if isfield(geom, 'iPer')
  geom.iPer = I(geom.iPer);
end

% YS: manily change here!!!!!!!!!!!!!!!
% geom.iBdry = I(geom.iBdry);

Idx = I(geom.iBdry)';
X = geom.X(:,Idx);
uIdx = find(X(2,:)>0 | abs(X(2,:))<1e-5);
lIdx = setdiff(1:length(Idx), uIdx);
p_up = X(1,uIdx);
[~, uid] = sort(p_up, 'descend' );
p_lw = X(1,lIdx);
[~, lid] = sort(p_lw, 'ascend' );
id = [uIdx(uid), lIdx(lid)];
geom.iBdry = Idx(id);



% 
% % if the geom structure has already been analyzed, rerun it
% if isfield(geom, 'geom_analyze')
%   if geom.geom_analyze
%     geom.geom_analyze = false;
%     geom = geom_analyze(geom);
%   end
% end


end


%% test routine
function test_geom_create_vacancies()

geom = geom_2dtri_longhex(2, 10, 10, 1.2, 'dir');
geom = geom_create_vacancies(geom, [1, 5, 6 7, 10]);
geom.plot(geom);

% geom = geom_2dtri_hexagon(30, 10, 1.2, 'per');
% geom = geom_create_vacancies(geom, [5, 10, 35]);
% geom.plot(geom);

% modify test_geom_analyze to check how it behaves after vacancies are
% removed

end
