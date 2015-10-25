function sos_region()

checkDependency('spotless');
checkDependency('mosek');

x = linspace(0, 2.5);
f = @(x) (x-1).^4 - (x-1).^2 + 1;

clf()
plot(x, f(x))
hold on

xx = linspace(0.1, 2.5, 15);
ok = false(size(xx));
for i = 1:length(xx)
  ok(i) = check_membership(xx(i));
  if ok(i)
    plot(xx(i), 0.9, 'go')
  else
    plot(xx(i), 0.9, 'rx')
  end

end

end

function ok = check_membership(x0)
  prog = spotsosprog;

  x = msspoly('x');
  prog = prog.withIndeterminate(x);

  v = monomials(x, 0:6);
  [prog, lambda1_coeffs] = prog.newFree(length(v));
  [prog, a] = prog.newFree(1);
  [prog, b] = prog.newFree(1);
  [prog, c] = prog.newFree(1);

  lambda1 = lambda1_coeffs'*v;

  V = (x - 1)^4 - (x - 1)^2 + 1;

  % P = (0.9 - V) + lambda1*(x - x0)^2 - 1e-4;
  P = (0.9 - V) + lambda1 - 1e-4;
  prog = prog.withSOS(P);
  prog = prog.withSOS(lambda1);
  prog = prog.withSOS(2*(x - x0) - lambda1 + (x0 - x));
  prog = prog.withSOS(2*(x0 - x) - lambda1 + (x - x0));
  % prog = prog.withSOS(a*x^2 + b*x + c);
  % prog = prog.withPos(a * x0 * x0 - 1);
  % prog = prog.withPos(2 * a * x0 + b);
  % prog = prog.withPos(a);
  prog = prog.withEqs(a*x0^2 + b*x0 + c);

  result = prog.minimize(0, @spot_mosek)
  ok = result.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE

  if ok
    result.eval(lambda1_coeffs)
  end

end

function [xstar, ok] = find_closest(x0)
  prog = spotsosprog;

  x = msspoly('x');
  prog = prog.withIndeterminate(x);

  v = monomials(x, 0:4);
  [prog, lambda1_coeffs] = prog.newFree(length(v));
  [prog, y] = prog.newFree(1);
  % [prog, a] = prog.newFree(1);

  lambda1 = lambda1_coeffs'*v;

  V = (x - 1)^4 - (x - 1)^2 + 1;

  P = lambda1*(0.9 - V) + (1/x0^2) * x^2 - 2*(1/x0^2)*x0*x + 1 - 1e-4;
  prog = prog.withSOS(P);
  prog = prog.withSOS(lambda1);
  % prog = prog.withEqs(a * x0 * x0 - 1);
  % prog = prog.withPos(a);

  result = prog.minimize(0, @spot_mosek)
  ok = result.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE

  result.eval(lambda1_coeffs)
end

