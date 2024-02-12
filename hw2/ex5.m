clc;
clear;
close all;

%% 5

%% 5.c Continuous time linearized model
disp('5.c');
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0];
u0 = [0; 0; 0];
syms J1 J2 J3 phi theta psi omega1 omega2 omega3 r1 r2 r3 M1 M2 M3;
x = [phi; theta; psi; omega1; omega2; omega3; r1; r2; r3];
u = [M1; M2; M3];
J = [J1; J2; J3];
f = dxdtwJ(x, u, J);
A = jacobian(f, x);
B = jacobian(f, u);
Ac = subs(A, [x; u], [x0; u0]);
Bc = subs(B, [x; u], [x0; u0]);

disp('Ac = ');
disp(Ac);
disp('Bc = ');
disp(Bc);

%% 5.d Discrete time linearized model
disp('5.d');
Ts = 2;
Ac = double(subs(A, [x; u; J], [x0; u0; [120; 100; 80]]));
Bc = double(subs(B, [x; u; J], [x0; u0; [120; 100; 80]]));
Cc = zeros(size(Ac, 1), size(Ac, 2));
Dc = zeros(size(Cc, 1), size(Bc, 2));
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, Ts, 'zoh');
% check if Ad is schur
disp('Eigenvalues of Ad = ');
disp(eig(Ad));

if all(abs(eig(Ad)) < 1)
  disp('Ad is schur');
else
  disp('Ad is not schur');
end

% check if Ad, Bd is controllable
% we dont care about controlling the r1, r2, r3 states, because they are reference
% hence ignore the last 3 states
if rank(ctrb(Ad, Bd)) == size(Ad, 1) - 3
  disp('(Ad, Bd) is controllable');
else
  disp('(Ad, Bd) is not controllable');
end

%%
% Ad is not a schur matrix since the eigenvalues are not inside the unit circle
%
% Ad, Bd is controllable since the rank of the controllability matrix is equal to the number of states that we want to control (6)

%% 5.e Implement terminal cost function
disp('5.e');
Q = diag([1, 1, 1, 0.01, 0.01, 0.01]);
R = diag([1, 1, 1]) * 0.01;
Ad = Ad(1:6, 1:6);
Bd = Bd(1:6, :);
global Pinf;
[Kinf, Pinf, e] = dlqr(Ad, Bd, Q, R);
disp('Pinf = ');
disp(Pinf);

%% 5.f Simulate the system
disp('5.f');
import casadi.*
mpc = import_mpctools();
Nx = 9;
Nu = 3;
Nt = 30;
Delta = 2;
Nsim = 100;
N = struct('x', Nx, 'u', Nu, 't', Nt);
x = NaN(Nx, Nsim + 1);
x(:, 1) = [-0.4; -0.8; 1.2; -0.02; -0.02; 0.02; 0; 0; 0];
u = NaN(Nu, Nsim);

for i = 0:1

  for t = 0:Nsim

    if i == 0
      odeFun = @(x, u) ode(x, u); % nonlinear
    else
      odeFun = @(x, u) (Ac * x + Bc * u); % linear
    end

    f = mpc.getCasadiFunc(odeFun, [Nx, Nu], {'x', 'u'}, 'rk4', true, 'Delta', Delta, 'M', 1);
    l = mpc.getCasadiFunc(@stagecost, [Nx, Nu], {'x', 'u'}, {'l'});
    Vf = mpc.getCasadiFunc(@termcost, [Nx], {'x'}, {'Vf'});
    lbx = -inf * ones(Nx, Nt + 1);
    lbx(4:6, :) = -0.03;
    ubx = inf * ones(Nx, Nt + 1);
    ubx(4:6, :) = 0.03;
    lbu = -0.1 * ones(Nu, Nt);
    ubu = 0.1 * ones(Nu, Nt);
    commonargs = struct('l', l, 'Vf', Vf, 'lb', struct('u', lbu, 'x', lbx), 'ub', struct('u', ubu, 'x', ubx));
    solvers = mpc.nmpc('f', f, 'N', N, 'Delta', Delta, '**', commonargs);

    if t < Nsim / 2
      r = [0; 0; 0];
    else
      r = [0.5; 0; -0.5];
    end

    x(7:9, t + 1) = r;
    solvers.fixvar('x', 1, x(:, t + 1));
    solvers.solve();
    % fprintf('%d: %s\n', t, solvers.status);
    if ~isequal(solvers.status, 'Solve_Succeeded')
      warning('%s failed at time %d', solvers.name, t);
    end

    solvers.saveguess();
    u(:, t + 1) = solvers.var.u(:, 1);
    x(:, t + 2) = solvers.var.x(:, 2);
  end

  % plot the figures
  set(0, 'DefaultLineLineWidth', 1.5);
  Time = 0:Delta:(Nsim + 1) * Delta;

  figure(i * 3 +1);

  if i == 0
    sgtitle("Euler Angles Nonlinear");
  else
    sgtitle("Euler Angles Linear");
  end

  subplot(3, 1, 1)
  plot(Time, x(1, :))
  grid on;
  hold on
  plot(Time, x(7, :), '-.')
  set(gca, 'xticklabel', [])
  ylabel('phi')
  xlim([Time(1), Time(end)]);
  legend('MPC', 'setpoint', 'Fontsize', 12)

  subplot(3, 1, 2)
  plot(Time, x(2, :))
  grid on;
  hold on
  plot(Time, x(8, :), '-.')
  set(gca, 'xticklabel', [])
  ylabel('theta')
  xlim([Time(1), Time(end)]);

  subplot(3, 1, 3)
  plot(Time, x(3, :))
  grid on;
  hold on
  plot(Time, x(9, :), '-.')
  set(gca, 'xticklabel', [])
  ylabel('psi')
  xlim([Time(1), Time(end)]);
  snapnow;

  figure(i * 3 + 2);

  if i == 0
    sgtitle("Angular velocities Nonlinear")
  else
    sgtitle("Angular velocities Linear")
  end

  subplot(3, 1, 1)
  plot(Time, x(4, :))
  grid on;
  hold on
  set(gca, 'xticklabel', [])
  ylabel('w1')
  xlim([Time(1), Time(end)]);
  yline(0.03, '--', 'LineWidth', 1);
  yline(-0.03, '--', 'LineWidth', 1);
  legend('MPC', 'Upper bound', "Lower bound", 'Fontsize', 12);

  subplot(3, 1, 2)
  plot(Time, x(5, :))
  grid on;
  hold on
  set(gca, 'xticklabel', [])
  ylabel('w2')
  xlim([Time(1), Time(end)]);
  yline(0.03, '--', 'LineWidth', 1);
  yline(-0.03, '--', 'LineWidth', 1)

  subplot(3, 1, 3)
  plot(Time, x(6, :))
  grid on;
  hold on
  set(gca, 'xticklabel', [])
  ylabel('w3')
  xlim([Time(1), Time(end)]);
  yline(0.03, '--', 'LineWidth', 1);
  yline(-0.03, '--', 'LineWidth', 1)
  snapnow;

  figure(i * 3 + 3);

  if i == 0
    sgtitle("Inputs Nonlinear")
  else
    sgtitle("Inputs Linear")
  end

  subplot(3, 1, 1)
  plot(Time(2:end), u(1, :));
  grid on;
  hold on
  set(gca, 'xticklabel', [])
  ylabel('M1')
  xlim([Time(1), Time(end)]);
  yline(0.1, '--', 'LineWidth', 1);
  yline(-0.1, '--', 'LineWidth', 1)
  legend('MPC', 'Upper bound', "Lower bound", 'Fontsize', 12)

  subplot(3, 1, 2)
  plot(Time(2:end), u(2, :))
  grid on;
  hold on;
  set(gca, 'xticklabel', []);
  ylabel('w2');
  xlim([Time(1), Time(end)]);
  yline(0.1, '--', 'LineWidth', 1);
  yline(-0.1, '--', 'LineWidth', 1)

  subplot(3, 1, 3)
  plot(Time(2:end), u(3, :))
  grid on;
  hold on
  set(gca, 'xticklabel', [])
  ylabel('w3')
  xlim([Time(1), Time(end)]);
  yline(0.1, '--', 'LineWidth', 1);
  yline(-0.1, '--', 'LineWidth', 1)
  snapnow;
end

%% 5.a Implement continuous time system function
function dxdt = ode(x, u)
  J = [120; 100; 80];
  dxdt = dxdtwJ(x, u, J);
end

function dxdtwJ = dxdtwJ(x, u, J)
  J1 = J(1);
  J2 = J(2);
  J3 = J(3);
  phi = x(1);
  theta = x(2);
  % psi = x(3);
  omega1 = x(4);
  omega2 = x(5);
  omega3 = x(6);
  % r1 = x(7);
  % r2 = x(8);
  % r3 = x(9);
  M1 = u(1);
  M2 = u(2);
  M3 = u(3);
  dstatedt = 1 / cos(theta) * ...
    [
   cos(theta) sin(phi) * sin(theta) cos(phi) * sin(theta);
   0 cos(phi) * cos(theta) -sin(phi) * cos(theta);
   0 sin(phi) cos(theta);
   ] * ...
    [omega1; omega2; omega3];
  dinputdt = [
              ((J2 - J3) * omega2 * omega3) / J1 + M1 / J1;
              ((J3 - J1) * omega3 * omega1) / J2 + M2 / J2;
              ((J1 - J2) * omega1 * omega2) / J3 + M3 / J3;
              ];
  dxdtwJ = [
            dstatedt;
            dinputdt;
            [0; 0; 0];
            ];
end

%% 5.b Implement stage cost function
function l = stagecost(x, u)
  r = x(7:9);
  Q = diag([1, 1, 1, 0.01, 0.01, 0.01]);
  R = diag([1, 1, 1]) * 0.01;
  e = [x(1:3) - r(:); x(4:6)];
  l = e' * Q * e + u' * R * u;
end

%% 5.e Implement terminal cost function
function Vf = termcost(x)
  r = x(7:9);
  e = [x(1:3) - r(:); x(4:6)];
  global Pinf;
  Vf = e' * Pinf * e;
end
