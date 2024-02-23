%% AE740 HW1 akshatdy

clc;
clear;
close all;

%% 2. Spacecraft attitude dynamics

%% 2.a Linearized dynamics

% differentiate to get jacobian
syms phi theta psi omega1 omega2 omega3 M1 M2 M3;
x = [phi; theta; psi; omega1; omega2; omega3; ];
u = [M1; M2; M3];
f = scdynamics(0, x, u);
A = jacobian(f, x);
B = jacobian(f, u);
% point of linearization
x0 = [0; 0; 0; 0; 0; 0; ];
nx = length(x0);
u0 = [0; 0; 0];
nu = length(u0);
Ac = double(subs(A, [x; u], [x0; u0]));
Bc = double(subs(B, [x; u], [x0; u0]));
% display
% disp('Ac = ');
% disp(Ac);
% disp('Bc = ');
% disp(Bc);
% convert to discrete time
Ts = 2;
Cc = zeros(size(Ac, 1), size(Ac, 2));
Dc = zeros(size(Cc, 1), size(Bc, 2));
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, Ts, 'zoh');
disp('Ad = ');
disp(Ad);
disp('Bd = ');
disp(Bd);

%% 2.b Simulate using MPT Toolbox

% MPT stuff
% tbxmanager restorepath
% mpt_init
% consts
x0 = [-0.4; -0.8; 1.2; -0.02; -0.02; 0.02];
Q = diag([1, 1, 1, 0.01, 0.01, 0.01]);
R = diag([0.01, 0.01, 0.01]);
% get P with DARE
[Kinf, Pinf, CLeig] = dlqr(Ad, Bd, Q, R);
% setup model
model = LTISystem('A', Ad, 'B', Bd, 'C', Cd, 'D', Dd, 'Ts', Ts);
model.u.min = -0.1 * ones(nu, 1);
model.u.max = 0.1 * ones(nu, 1);
model.u.penalty = QuadFunction(R);
model.x.penalty = QuadFunction(Q);
model.x.with('terminalPenalty');
model.x.terminalPenalty = QuadFunction(Pinf);

Tsim = 200;
times = 0:Ts:Tsim;
h = Ts / 10;

ctrl2 = MPCController(model, 2);
ctrl20 = MPCController(model, 20);
data.X2 = zeros(nx, Tsim / Ts);
data.X20 = zeros(nx, Tsim / Ts);
data.U2 = zeros(nu, Tsim / Ts);
data.U20 = zeros(nu, Tsim / Ts);
data.execTime2 = zeros(1, Tsim / Ts);
data.execTime20 = zeros(1, Tsim / Ts);
x2 = x0;
x20 = x0;
dataIdx = 1;

for t = times
  tic;
  u2 = ctrl2.evaluate(x2);
  data.execTime2(dataIdx) = toc;
  data.U2(:, dataIdx) = u2;
  [~, Xode2] = ode45(@ (t, x) scdynamics(t, x, u2), [t + h:h:t + Ts], x2);
  x2 = Xode2(end, :)';
  data.X2(:, dataIdx) = x2;

  tic;
  u20 = ctrl20.evaluate(x20);
  data.execTime20(dataIdx) = toc;
  data.U20(:, dataIdx) = u20;
  [~, Xode20] = ode45(@ (t, x) scdynamics(t, x, u20), [t + h:h:t + Ts], x20);
  x20 = Xode20(end, :)';
  data.X20(:, dataIdx) = x20;

  dataIdx = dataIdx + 1;
end

figure();
sgtitle("MPT");

for i = 1:3
  subplot(3, 3, i)
  plot(times, data.X2(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.X20(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("x_" + num2str(i));
  legend show
end

for i = 4:6
  subplot(3, 3, i)
  plot(times, data.X2(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.X20(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("w_" + num2str(i));
  legend show
end

for i = 1:3
  subplot(3, 3, 6 + i)
  plot(times, data.U2(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.U20(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("u_" + num2str(i));
  legend show
end

figure();
sgtitle('Time to compute control input');
plot(times, data.execTime2, 'DisplayName', 'N=2');
grid on;
hold on;
plot(times, data.execTime20, 'DisplayName', 'N=20');
xlabel('t [s]');
ylabel('Time [s]');
legend show

function [xdot] = scdynamics(t, x, u)
  J1 = 120;
  J2 = 100;
  J3 = 80;
  phi = x(1);
  theta = x(2);
  % psi = x(3);
  omega1 = x(4);
  omega2 = x(5);
  omega3 = x(6);
  M1 = u(1);
  M2 = u(2);
  M3 = u(3);
  angledot = 1 / cos(theta) * ...
    [
   cos(theta) sin(phi) * sin(theta) cos(phi) * sin(theta);
   0 cos(phi) * cos(theta) -sin(phi) * cos(theta);
   0 sin(phi) cos(phi);
   ] * ...
    [omega1; omega2; omega3];
  velocitydot = [
                 ((J2 - J3) * omega2 * omega3) / J1 + M1 / J1;
                 ((J3 - J1) * omega3 * omega1) / J2 + M2 / J2;
                 ((J1 - J2) * omega1 * omega2) / J3 + M3 / J3;
                 ];
  xdot = [
          angledot;
          velocitydot;
          ];
end
