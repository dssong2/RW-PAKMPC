function draw_combined(ax, q, p)
%DRAW_COMBINED  Draw the combined spacecraft of build_pakmpc_multibody at attitude q (scalar first).
%   [guess] Same geometry as build_pakmpc_multibody: servicer cube, target box docked on the
%   servicer +y face through a ring, two three-link arms and two target solar panels. Body points
%   are taken about the combined center of mass and rotated to the inertial frame by R(q).

ls = sqrt(6*p.js(1, 1)/p.ms);
lt = sqrt(6/p.mt*(trace(p.jt) - 2*diag(p.jt)'));
r = [0; ls/2 + 0.2 + lt(2)/2; 0];
c = p.mt*r/(p.ms + p.mt);
rq = q2r(q);
hold(ax, 'on');
box3(ax, rq, -c, [ls ls ls], [0.75 0.75 0.8]);                    % servicer
box3(ax, rq, r - c, lt, [0.85 0.7 0.3]);                          % target
for s = [-1 1]
    box3(ax, rq, r - c + [s*(lt(1)/2 + 0.25 + 0.9); 0; 0], [1.8 0.9 0.03], [0.1 0.15 0.45]);
    pts = [s*0.3*ls, ls/2, 0; s*(lt(1)/2 + 0.35), ls/2 + 0.25, 0;
           s*(lt(1)/2 + 0.45), r(2) - 0.2, 0; s*lt(1)/2, r(2), 0]' - c;
    w = rq*pts;
    plot3(ax, w(1, :), w(2, :), w(3, :), '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 3);
    b = rq*([s*lt(1)/2 s*(lt(1)/2 + 0.25); 0 0; 0 0] + r - c);
    plot3(ax, b(1, :), b(2, :), b(3, :), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
end
g = rq*([0 0; ls/2 ls/2 + 0.2; 0 0] - c);
plot3(ax, g(1, :), g(2, :), g(3, :), '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 8);
% body axes
for i = 1:3
    e = rq(:, i)*3.2;
    quiver3(ax, 0, 0, 0, e(1), e(2), e(3), 0, 'Color', [0.8 0 0]*(i == 1) + [0 0.6 0]*(i == 2) + ...
        [0 0 0.9]*(i == 3), 'LineWidth', 1.5, 'MaxHeadSize', 0.4);
end
end

function box3(ax, rq, ctr, dims, col)
[x, y, z] = ndgrid([-1 1]/2, [-1 1]/2, [-1 1]/2);
v = rq*([x(:)*dims(1) y(:)*dims(2) z(:)*dims(3)]' + ctr);
f = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
patch(ax, 'Vertices', v', 'Faces', f, 'FaceColor', col, 'EdgeColor', 0.3*col, 'FaceLighting', 'gouraud');
end
