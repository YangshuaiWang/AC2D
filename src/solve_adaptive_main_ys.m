
function [geom, model, Y, adapflow] = solve_adaptive_main_ys(geom, model, Y, dirname, sample, fn, adapflow)
% check input
N = geom.N;
dof = geom.nX;

% current adaptive stage
adapk = geom.adapk;

% refine number
refk = 3;

load(['./out/', dirname, '/setup'], 'data');
F0 = data.F0;

Sfx = sprintf('%d-%d.mat', N, dof);

str = [fn, '/', Sfx];
% load  or compute related data.
slflag = 0;
if exist(str, 'file')
    load(str);
    slflag = 1;
end

geom = geom_2dtri_vacTri(geom);

if slflag == 1 && exist('estimate', 'var')
    fprintf('Data: %s loaded!\n', str);
else
%   [estimate, ~] = solve_adaptive_error(geom, model, Y, dirname, F0);
%   with a given stability constant, call solve_adaptive_error as
%   below! Or above to compute stab structure.
    geom = geom_2dtri_h(geom);
    [~, ik] = sort(geom.volT, 'descend');
    iT_refine = ik(1:refk);
    xh = geom.heart(:,iT_refine);
    adapflow(adapk).refX = xh;
    estimate.T_refine = iT_refine; % = solve_adaptive_error(geom, model, Y, fn); # only need to change here!
end

opt = 'on';
fprintf('Refining selected elements...\n')

disp('Element refinement...');

iT_rfn = estimate.T_refine;

%%%%%%% construct new geom %%%%%%%%%%%%%%
geom_new = geom_2dtri_mcrack(model.Kc, geom.K+1, geom.N, model.rCutH, geom.alpha);
% geom_new = geom_analyze_all(geom_new);
for j = 1:adapk
    xhj = adapflow(j).refX;
    iT_rfn = tsearchn(geom_new.X', geom_new.T', xhj');
    [row, ~] = size(iT_rfn);
    if row == 1; iT_rfn = iT_rfn'; end;
    geom_new = solve_adaptive_refine(geom_new, opt, sample, iT_rfn);
%     geom_new = geom_analyze_all(geom_new);
end
geom = geom_analyze_all(geom_new);
% geom = geom_new;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% iX = find(geom.di==0);
% iAtom = find(geom.volX ~= 0);
% iT_atom = geom_2dtri_aIXT(geom, iAtom);
% iT_intf = geom_2dtri_aIXT(geom, iX);
% iT_exd = union(iT_atom, iT_intf);
% iT = intersect(iT_rfn, iT_exd);
% % proportion = sum(estimate.err(iT))/estimate.total;
% % iT_rfn = setdiff(iT_rfn, iT_atom);
% iT_rfn = setdiff(iT_rfn, iT_exd)';
% % if proportion >= 0.05;
% %     reply = 'Y';
% % else
% %     reply = 'N';
% % end
% 
% fprintf(['-----------------------------------------------------' ...
%              '-----------------------\n']);
% fprintf('\t#\testimator = %3.5f\n', estimate.total);
% fprintf(fid, '\t#\testimator = %3.5f\n', estimate.total);
% fclose(fid);
% 
% if numel(iT_rfn) == 1; iT_rfn = [iT_rfn; iT_rfn]; end
% % if isempty(iT_rfn) && reply == 'N'
% %     disp('No element to refine');
% %     return;
% % end
% % ======================================================================= %
% % TODO: not right now
% % if strcmpi(reply, 'Y')
% %     sample = lower(sample);
% %     switch sample
% %         case 'zze'
% %             [geom_ex, iT_ex] = solve_adaptive_ExAtom(geom, 'geom-ex', error);
% %         case 'amc'
% %             [geom_ex, iT_ex] = solve_adaptive_ExAtom(geom, 'geom-ex', error);
% % %             geom = geom_2dtri_jigglemesh(geom);
% %             geom_ex = gqc23_prep_geom(geom);
% %         case 'asv'
% %             [geom_ex, iT_ex] = solve_adaptive_ExAtom(geom, 'geom-ex', error);
% % %             geom = geom_2dtri_jigglemesh(geom);
% %             geom_ex = gqc23_prep_geom(geom_ex);
% %         case 'amv'
% %             geom_ex = solve_multiVc_ExAtom(geom);
% %         otherwise
% %             warning('unkown sampe type! treat as Amcrack');
% %             geom_ex = solve_adaptive_ExAtom(geom, 'geom-ex');
% %     end
% % end
% % ======================================================================= %
% % figure(1);
% % clf
% % error_distribute(geom, estimate.err, jet)
% 
% estimate.T_refine = iT_rfn;
% estimate.total = sum(estimate.err(iT_rfn));
% 
% % save(str, 'estimate', '-append')
% % if strcmpi(reply, 'Y')
% %     [geom_ex, iT_ex] = solve_adaptive_ExAtom(geom, 'geom-ex', estimate, fn);
% %     geom_ex = gqc23_prep_geom(geom_ex);
% % else 
% %     geom_ex = geom;
% %     iT_ex = [];
% % end
% 
% [geom_ex, iT_ex] = solve_adaptive_ExAtom(geom, 'geom-ex', estimate, fn);
% geom_ex = gqc23_prep_geom(geom_ex);
% 
% geom = geom_2dtri_h(geom);
% geom_ex = geom_2dtri_h(geom_ex);
% heart = geom.heart;
% iT_rfn = setdiff(iT_rfn, iT_ex);
% K1 = geom.K1;
% % geom = geom_ex;
% geom_ex = geom_2dtri_vacTri(geom_ex);
% heart = heart(:,iT_rfn);
% iT_rfn = tsearchn(geom_ex.X', geom_ex.Tri, heart');
% % iT_rfn = geom_2dtri_aPT(heart, geom_ex);
% iT_rfn = iT_rfn';
% geom = geom_ex;
% iT_small = find(geom.volT<0.5+1e-5)';
% iT_rfn = setdiff(iT_rfn, iT_small);
% 
% if isempty(iT_rfn) && isempty(iT_ex)
%     disp('No element to refine');
%     return;
% end
% % iX = find(geom.di==0);
% % iAtom = find(geom.volX ~= 0);
% % iT_atom = geom_2dtri_aIXT(geom, iAtom);
% % iT_intf = geom_2dtri_aIXT(geom, iX);
% % iT_exd = union(iT_atom, iT_intf);
% % iT_rfn = setdiff(iT_rfn, iT_exd)';
% %% feature to-be refined elements adjacent to interface,
% iX = find(geom.di == 0);
% iT_old = geom_2dtri_aIXT(geom, iX);
% % iAtom = find(geom.volX ~= 0);
% % iT_atom = geom_2dtri_aIXT(geom, iAtom);
% 
% E = geom.aTE(:, iT_old);
% E = unique(E);
% iT_new = geom.aET(:, E);
% iT_new = unique(iT_new);
% % iT_new = setdiff(iT_new, iT_atom);
% % iT_new = intersect(iT_new, iT_rfn);
% 
% iT_rfn = setdiff(iT_rfn, iT_new);
% 
% fid = fopen([fn, '/progress.log'], 'a');
% fprintf('\t#\t#(Bisect) = %d\n', numel(iT_rfn));
% fprintf(fid, '\t#\t#(Bisect) = %d\n', numel(iT_rfn));
% fclose(fid);
% 
% if numel(iT_rfn) == 1; iT_rfn = [iT_rfn; iT_rfn]; end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isempty(iT_rfn)
%     fprintf('No element to refine!\n');
%     return;        
% else
%     [row, ~] = size(iT_rfn);
%     if row == 1; iT_rfn = iT_rfn'; end;
%     geom = solve_adaptive_refine(geom, opt, sample, iT_rfn);
%     geom = geom_analyze(geom, 1);
%     geom = gqc23_prep_geom(geom);
%     geom = geom_2dtri_h(geom);
% end

figure(adapk)
geom.plot(geom)
% Sf = sprintf('%d.fig', adapk);
% str = [fn, '/geom-', Sf];
% savefig(str)

next = geom.nX;
save(str, 'next');

N = geom.N;
dof = geom.nX;
Sfx = sprintf('%d-%d.mat', N, dof);
% % sfx = sprintf('%d_%d', N, dof);

str = [fn, '/soln-', Sfx];

if exist(str, 'file')
    load(str, 'geom', 'Y', 'E');
    fprintf('Data: %s loaded!\n', str);
else
    [Y, E, geom] = solve_adaptive_QcSoln(geom, model, F0);
    save(str, 'geom', 'Y', 'E', 'model');
end
return
end
