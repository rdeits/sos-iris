function sos_sampling()

checkDependency('spotless');
checkDependency('mosek');


xy_frame = CoordinateFrame('xy', 2);
xy = xy_frame.getPoly();

true_obstacles = {xy(1) - 2,...
             -xy(1) - 2,...
             xy(2) - 2,...
             -xy(2) - 2,...
             -xy(1)^2 - xy(2)^2 + 1,...
             -2*(xy(1)+1.5)^2 - 2*xy(2)^2 + 1};
known_obstacles = {xy(1) - 2,...
             -xy(1) - 2,...
             xy(2) - 2,...
             -xy(2) - 2};

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
% [prog, Obs_coeffs] = prog.newFree(length(mon));
% Obs = Obs_coeffs' * mon;

% sample_x = linspace(-2, 2, 17);
% sample_y = linspace(-2, 2, 17);
% [sample_x, sample_y] = meshgrid(sample_x, sample_y);
% sample_points = [sample_x(:)'; sample_y(:)'];
num_samples = 200;
sample_points = zeros(2, num_samples);
for j = 1:num_samples
  sample_points(:,j) = random('uniform', [-2.5; -2.5], [2.5; 2.5]);
end
% for j = 1:length(obstacles)
%   sample_points = sample_points(:, msubs(obstacles{j}, xy, sample_points) <= 0);
% end


[prog, sample_costs] = prog.newFree(size(sample_points, 2));

for j = 1:size(sample_points, 2)
  known_to_be_in_obstacle = false;
  in_obstacle = false;
  for k = 1:length(known_obstacles)
    if (msubs(known_obstacles{k}, xy, sample_points(:,j)) > 0)
      known_to_be_in_obstacle = true;
      break;
    end
  end
  for k = 1:length(true_obstacles)
    if msubs(true_obstacles{k}, xy, sample_points(:,j)) > 0
      in_obstacle = true;
      break;
    end
  end
  if in_obstacle && ~known_to_be_in_obstacle
    prog = prog.withPos(msubs(V, xy, sample_points(:,j)));
  end
  % elseif ~in_obstacle
  prog = prog.withPos(sample_costs(j)-msubs(V, xy, sample_points(:,j)));
  % end
  prog = prog.withPos(sample_costs(j));
end

for j = 1:length(known_obstacles)
  prog = prog.withSOS(V - known_obstacles{j});
end
% prog = prog.withSOS(V - Obs);

t0 = tic();
result = prog.minimize(sum(sample_costs), @spot_mosek)
toc(t0)
V = result.eval(V)
% Obj = result.eval(Obs)

figure(1);
clf
hold on

for j = 1:size(sample_points, 2)
  if msubs(V, xy, sample_points(:,j)) <= 0
    color = 'g';
  else
    color = 'k';
    for k = 1:length(true_obstacles)
      if msubs(true_obstacles{k}, xy, sample_points(:,j)) >= 0
        color = 'r';
        break;
      end
    end
  end
  plot(sample_points(1,j), sample_points(2,j), 'o', 'Color', color, 'MarkerSize', 10, 'MarkerFaceColor', color);
end

[X, Y] = meshgrid(linspace(-2.2, 2.2), linspace(-2.2, 2.2));
Z = reshape(msubs(V, xy, [X(:)'; Y(:)']), size(X));
contour(X, Y, Z, [0, 0]);



figure(2);
clf
hold on
% C = 2*ones(size(X));
% C(Z <= 0) = 1;
surf(X, Y, Z);
colormap hsv

% for j = 1:length(obstacles)
%   Z = reshape(msubs(obstacles{j}, xy, [X(:)'; Y(:)']), size(X));
%   C = -0.5*ones(size(X));
%   C(Z >= 0) = -2;
%   mesh(X, Y, Z, C);
% end
% Z = reshape(msubs(Obs, xy, [X(:)'; Y(:)']), size(X));
% C = -0.5*ones(size(X));
% C(Z >= 0) = -2;
% mesh(X, Y, Z, C);

% colormap hsv
% zlim([-0.5, 2.5]);

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
