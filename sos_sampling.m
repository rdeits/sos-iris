function sos_sampling()

checkDependency('spotless');
checkDependency('mosek');

% Our free space is the set [-2, -1] + [1, 2]
x = msspoly('x', 1);
obstacles = {x - 2,...
             -x - 2,...
             -x^2 + 1};

prog = spotsosprog();
prog = prog.withIndeterminate(x);
[prog, V_coeffs] = prog.newFree(5);
mon = monomials(x, 0:4);
V = V_coeffs' * mon;

sample_points = linspace(-2, 2, 10);
[prog, sample_costs] = prog.newFree(numel(sample_points));

for j = 1:length(sample_points)
  prog = prog.withPos(sample_costs(j)-msubs(V, x, sample_points(j)));
  prog = prog.withPos(sample_costs(j));
end

for j = 1:length(obstacles)
  prog = prog.withSOS(V - obstacles{j});
end

result = prog.minimize(sum(sample_costs), @spot_mosek)
V = result.eval(V)

figure(1);
clf
hold on
plot(sample_points, zeros(size(sample_points)), 'kx');
for j = 1:length(obstacles)
  xx = linspace(-3, 3);
  yy = msubs(obstacles{j}, x, xx);
  plot(xx, yy, 'r-');
end
xx = linspace(-3, 3);
yy = msubs(V, x, xx);
plot(xx, yy, 'g-');

end
