function sos_sampling()

checkDependency('spotless');
checkDependency('mosek');


xy_frame = CoordinateFrame('xy', 2);
xy = xy_frame.getPoly();

obstacles = {xy(1) - 2,...
             -xy(1) - 2,...
             xy(2) - 2,...
             -xy(2) - 2,...
             -xy(1)^2 - xy(2)^2 + 1,...
             -2*(xy(1)+1.5)^2 - 2*xy(2)^2 + 1};

% x = msspoly('x', 1);
% y = msspoly('y', 1);
% obstacles = {x - 2,...
%              -x - 2,...
%              y - 2,...
%              -y - 2,...
%              -x^2 - y^2 + 1};

prog = spotsosprog();
prog = prog.withIndeterminate(xy);
V_degree = 6;
mon = monomials(xy, 0:V_degree);
[prog, V_coeffs] = prog.newFree(length(mon));
V = V_coeffs' * mon;

sample_x = linspace(-2, 2, 17);
sample_y = linspace(-2, 2, 17);
[sample_x, sample_y] = meshgrid(sample_x, sample_y);
sample_points = [sample_x(:)'; sample_y(:)'];
[prog, sample_costs] = prog.newFree(size(sample_points, 2));

for j = 1:size(sample_points, 2)
  prog = prog.withPos(sample_costs(j)-msubs(V, xy, sample_points(:,j)));
  prog = prog.withPos(sample_costs(j));
end

for j = 1:length(obstacles)
  prog = prog.withSOS(V - obstacles{j});
end

t0 = tic();
result = prog.minimize(sum(sample_costs), @spot_mosek)
toc(t0)
V = result.eval(V)
V_plus_1 = SpotPolynomialLyapunovFunction(xy_frame, V+1);

figure(1);
clf
hold on

% V_plus_1.plotFunnel(struct('x0', [1.5; 1.5], 'tol', 1e-4));

for j = 1:size(sample_points, 2)
  if msubs(V, xy, sample_points(:,j)) <= 0
    color = 'g';
  else
    color = 'k';
    for k = 1:length(obstacles)
      if msubs(obstacles{k}, xy, sample_points(:,j)) >= 0
        color = 'r';
        break;
      end
    end
  end
  plot(sample_points(1,j), sample_points(2,j), 'o', 'Color', color, 'MarkerSize', 10, 'MarkerFaceColor', color);
end

[X, Y] = meshgrid(linspace(-2.2, 2.2), linspace(-2.2, 2.2));
Z = reshape(msubs(V, xy, [X(:)'; Y(:)']), size(X));
figure(2);
clf
hold on
C = 0.5*ones(size(X));
C(Z <= 0) = 2;
mesh(X, Y, Z, C);

for j = 1:length(obstacles)
  Z = reshape(msubs(obstacles{j}, xy, [X(:)'; Y(:)']), size(X));
  C = -0.5*ones(size(X));
  C(Z >= 0) = -2;
  mesh(X, Y, Z, C);
end

colormap hsv
zlim([-0.5, 2.5]);

% plot(sample_points, zeros(size(sample_points)), 'kx');
% for j = 1:length(obstacles)
%   xx = linspace(-3, 3);
%   yy = msubs(obstacles{j}, x, xx);
%   plot(xx, yy, 'r-');
% end
% xx = linspace(-3, 3);
% yy = msubs(V, x, xx);
% plot(xx, yy, 'g-');

end
