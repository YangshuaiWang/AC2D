
function pub_vis_bqce(geom, Y, mrka, mrkc, lwa, lwc)

% split the triangulation into atomistic and continuum
iTa = zeros(1, geom.nT);
for k = 1:geom.nT
  if min( geom.volX(geom.T(:,k)) ) == 1
    iTa(k) = 1;
  end
end

Ta = geom.T(:, iTa == 1);
Tc = geom.T(:, iTa == 0);

if lwc > 0
  triplot(Tc', Y(1,:)', Y(2,:)', 'k', 'Linewidth', lwc); hold on;
end
if lwa > 0
  triplot(Ta', Y(1,:)', Y(2,:)', 'r', 'Linewidth', lwa); hold on;
end
if mrkc > 0
  plot(Y(1,geom.onL==false), Y(2,geom.onL==false), 'b.', ...
    'Markersize', mrkc); hold on;
end
if mrka > 0
  scatter(Y(1,:), Y(2,:), mrka*geom.volX+eps, geom.volX, 'filled');
%   plot(Y(1,geom.onL==true), Y(2,geom.onL==true), 'r.', ...
%     'Markersize', mrka); hold on;
end

axis equal;
hold off;

end


