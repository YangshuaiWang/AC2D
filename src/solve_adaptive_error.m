% function [estimate, stab] = solve_adaptive_error(geom, model, Y, dirname, F0)
% Program: For inputs, compute error, stability variables and display info
%       on screen.
%   Input: geom, proper geometry structure; model, proper energy structure;
%       Y, displacements.
%   Output: square, edges triangle structure; error, error structure; stab,
%       storage of stability variables; nrm, w1p norm of approxes to exact.
% Author: M.Liao
% Versiion: First release   11-Mar-2015
%           1.1 Applied the div-free tensor field to correct STh.
%                                                   Mar-10-2016
function [estimate, stab] = solve_adaptive_error(geom, model, Y, fn)
% % TODO, check inputs and test routine.
N = geom.N;
disp('Computing residual...');
% compute square
dof = size(Y,2);
Sfx = sprintf('%d-%d.mat', N, dof);

str = [fn, '/', Sfx];
% load  or compute related data.
existflag = 0;
if exist(str, 'file')
	existflag = 1;
	load(str);	
end

if existflag == 0 || exist('square', 'var') ==0
    fprintf(['-----------------------------------------------------' ...
                 '-----------------------\n']);
    fprintf('Computing square...\n');
	tic; 
	square = geom_2dtri_scross(geom);
    fprintf('Interpolating square...\n');
	square = geom_interp_scross(geom, square, Y, model);
    fprintf('Analyzing square...\n');
	square = geom_scross_neigs(geom, square, Y, model);
	toc;
	square.exT = toc;
    if existflag == 0;
        save(str, 'square');
        existflag = 1;
    else    
        save(str, 'square', '-append' );
    end
end
% compute estimate
if existflag == 0 || exist('estimate', 'var') ==0
    fprintf(['-----------------------------------------------------' ...
                 '-----------------------\n']);
	tic;
    fprintf('Computing stress tensor...\n');
    STh = get_stress_tensor(geom, model, Y);
    [~, STe] = scross_stress_tensor(square, model, geom, Y);
    [S, ufs, divFC] = error_div_free(geom, square, STh, STe);
    STh = gqc23_stress_tensor(STh, S, ufs, divFC);
    fprintf('Computing H1 estimator...\n');
	estimate = error_main(geom, square, STh, STe);
    etaT = estimate_truncation(geom, STh, model);
    estimate.truncat = etaT;
	toc;
	estimate.exT = toc;
    save(str, 'estimate', 'STh', 'STe', '-append');
end

fprintf('Total residual is %3.5f.\n', estimate.total);
if nargout > 1
    % compute stability variables
    fprintf(['-----------------------------------------------------' ...
                 '-----------------------\n']);
    if existflag == 0 || exist('stab', 'var') ==0
        tic;
        stab = stability_constant(geom, model, square, Y, dirname);
        toc;
        stab.exT = toc;
        save(str, 'stab', '-append' );
    end
    % This load part is temporary, may be regarded as augments in the future.
    estimate.Est = estimate.total/stab.const;
    save(str, 'estimate', '-append');
end
fprintf('   # Dof               : %d \n', geom.nX);
fprintf('   # estimator         : %3.5f \n', estimate.total);
fprintf('   # elapsed time      : %3.5f \n', square.exT+estimate.exT);
% fprintf('   # elapsed time      : %3.5f \n', square.exT+estimate.exT+stab.exT);
fprintf(['-----------------------------------------------------' ...
       '-----------------------\n']);
return
end
