% function [geom, existflag] = solve_adaptive_refine(geom, opt, sample, iT_rfn)
% Program: adapt geom according to K, T_rfn and opt, then compute QC
%       solution.
%   Input: geom, geometry struct; model, energy struct; K, number of
%       element layers to construct; T_rfn, indexes of elements to bisect; 
%       opt, option of T_rfn i.e. jiggle on or off.
%   Output: new geometry struct geom, unchanged model, and QC solution Y;
% Author: M.Liao
% Version: First release    12-Mar-2015
%       1.01 remove input K to do element refinement only. Mar-31-2015
%       useless T.T Mar-31-2015
%       1.02 move counting iT_rfn to solve_adaptive_main, hence reduce
%       input arguments to goem, opt, sample.       Mar-12-2016
function geom = solve_adaptive_refine(geom, opt, sample, iT_rfn)
% TODO: check input, test routine.
% mesh refinement info
% existflag = 1;
% disp('Element refinement...');
% T_rfn = error.T_refine;
% I_node = find(geom.di==0);
% T_intf = zeros(1,0);
% for i = I_node
%     t = ceil(find(geom.T(:)==i)/3);
%     t = t(geom.volT(t)>0.5);
%     T_intf = [T_intf, t'];
% end
% T_intf = [T_intf, square.iTI];
% T_intf = unique(T_intf);
% T_rfn = setdiff(T_rfn, T_intf)';
% if length(T_rfn) == 1; T_rfn = [T_rfn; T_rfn]; end
% if isempty(T_rfn);
%     disp('No element to refine');
%     existflag = 0;
%     return;
% end

fprintf(['-----------------------------------------------------' ...
       '-----------------------\n']);
% fprintf('     Refine elements indexes:\n');
fprintf('Bisecting elements...\n');
% fprintf('     ');
% for i = 1:length(iT_rfn)
%     fprintf('     %5d\t', iT_rfn(i));
%     if mod(i,6) == 0;
%         fprintf('\n');
%     end
% end
% fprintf('\n');
% fprintf(['-----------------------------------------------------' ...
%        '-----------------------\n']);

% Default value
Value = 'on';

if nargin >= 3
    Value = eval('opt');
    if ~strcmpi(Value,'off') && ~strcmpi(Value,'on');
            error(message('geom:refinemesh:OptInvalidString'));
    end  
end

% sample = lower(sample);
% switch sample
%     case 'zze'
%         geom = geom_2dtri_refinemesh(geom, iT_rfn, 'jiggle', Value);
%     case 'amc'
%         geom = geom_2dtri_refinemesh(geom, iT_rfn, 'jiggle', Value);
%     case 'asv'
%         geom = geom_2dtri_refinemesh(geom, iT_rfn, 'jiggle', Value);
%     case 'amv'
%         geom = geom_multiVc_refinemesh(geom, iT_rfn, 'jiggle', Value);
%     otherwise
%         warning('unkown sampe type! treat as Amcrack');
%         geom = geom_2dtri_refinemesh(geom, iT_rfn, 'jiggle', Value);
% end
geom = geom_2dtri_refinemesh(geom, iT_rfn, 'jiggle', Value);
% geom = geom_multiVc_refinemesh(geom, T_rfn, 'jiggle', Value);
return
end
