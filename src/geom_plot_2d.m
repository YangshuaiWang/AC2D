%% function geom_plot_2d(geom)
%
% plot the geometry and some geometry information
% TODO: cleanup, include more information to plot
%


% COMMENT: use geom_analyze to test this routine

function geom_plot_2d(geom)

triplot(geom.T', geom.X(1,:)', geom.X(2,:)'); hold on;
plot(geom.X(1,geom.onL==false), geom.X(2,geom.onL==false), 'b.', ...
  'Markersize', 10);
plot(geom.X(1,geom.onL==true), geom.X(2,geom.onL==true), 'r.', ...
  'Markersize', 10);
if isfield(geom, 'iBdry')
  plot(geom.X(1, geom.iBdry), geom.X(2, geom.iBdry), 'k*');
end
axis equal;
hold off;

end


