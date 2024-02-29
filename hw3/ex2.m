%% AE740 HW1 akshatdy

clc;
clear;
close all;
clear variables
format shortG;

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
scd = c2d(ss(Ac, Bc, Cc, Dc), Ts, 'zoh');
Ad = scd.A;
disp('Ad = ');
disp(Ad);
Bd = scd.B;
disp('Bd = ');
disp(Bd);
Cd = scd.C;
Dd = scd.D;

%% 2.b Simulate using MPT Toolbox
disp("2.b Simulate using MPT Toolbox");
% consts
x0 = [-0.4; -0.8; 1.2; -0.02; -0.02; 0.02];
Q = diag([1, 1, 1, 0.01, 0.01, 0.01]);
R = diag([0.01, 0.01, 0.01]);
% get P with DARE
[Kinf, Pinf, CLeig] = dlqr(Ad, Bd, Q, R);

% MPT stuff, do once
tbxmanager restorepath
mpt_init

% setup model
modelMPT = LTISystem('A', Ad, 'B', Bd, 'C', Cd, 'D', Dd, 'Ts', Ts);
modelMPT.u.min = -0.1 * ones(nu, 1);
modelMPT.u.max = 0.1 * ones(nu, 1);
modelMPT.u.penalty = QuadFunction(R);
modelMPT.x.penalty = QuadFunction(Q);
modelMPT.x.with('terminalPenalty');
modelMPT.x.terminalPenalty = QuadFunction(Pinf);
ctrl2 = MPCController(modelMPT, 2);
ctrl20 = MPCController(modelMPT, 20);

%% 2.c Simulate using Hybrid MPC Toolbox
disp("2.c Simulate using Hybrid MPC Toolbox");
%%
% Both the toolboxes give the same results, so we can use either one for
% simulation, and it assures us that the simulation is correct.

% Hybrid toolbox stuff
addpath( ...
  genpath('/home/akshatd/mine/umich/sem2/740/toolboxes/hybtbx-linux/'), ...
'-end')
clear cost constraints horizon2 horizon20;

% setup model
cost.Q = Q;
cost.R = R;
cost.P = Pinf;
constraints.umin = -0.1 * ones(nu, 1);
constraints.umax = 0.1 * ones(nu, 1);
NHoriz = 2;
horizon2.Nu = NHoriz;
horizon2.N = NHoriz;
horizon2.Ncu = NHoriz - 1;
horizon2.Ncy = NHoriz - 1;
ctrl2h = lincon(scd, 'reg', cost, horizon2, constraints, 'qpact', 0);

NHoriz = 20;
horizon20.Nu = NHoriz;
horizon20.N = NHoriz;
horizon20.Ncu = NHoriz - 1;
horizon20.Ncy = NHoriz - 1;
ctrl20h = lincon(scd, 'reg', cost, horizon20, constraints, 'qpact', 0);
% do sim for both toolboxes
Tsim = 200;
times = 0:Ts:Tsim;
h = Ts / 10;

data.X2 = zeros(nx, Tsim / Ts);
data.X20 = zeros(nx, Tsim / Ts);
data.X2h = zeros(nx, Tsim / Ts);
data.X20h = zeros(nx, Tsim / Ts);
data.U2 = zeros(nu, Tsim / Ts);
data.U20 = zeros(nu, Tsim / Ts);
data.U2h = zeros(nu, Tsim / Ts);
data.U20h = zeros(nu, Tsim / Ts);
data.execTime2 = zeros(1, Tsim / Ts);
data.execTime20 = zeros(1, Tsim / Ts);
data.execTime2h = zeros(1, Tsim / Ts);
data.execTime20h = zeros(1, Tsim / Ts);
x2 = x0;
x20 = x0;
x2h = x0;
x20h = x0;
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

  tic;
  u2h = eval(ctrl2h, x2h);
  data.execTime2h(dataIdx) = toc;
  data.U2h(:, dataIdx) = u2h;
  [~, Xode2h] = ode45(@ (t, x) scdynamics(t, x, u2h), [t + h:h:t + Ts], x2h);
  x2h = Xode2h(end, :)';
  data.X2h(:, dataIdx) = x2h;

  tic;
  u20h = eval(ctrl20h, x20h);
  data.execTime20h(dataIdx) = toc;
  data.U20h(:, dataIdx) = u20h;
  [~, Xode20h] = ode45(@ (t, x) scdynamics(t, x, u20h), [t + h:h:t + Ts], x20h);
  x20h = Xode20h(end, :)';
  data.X20h(:, dataIdx) = x20h;

  dataIdx = dataIdx + 1;
end

% some other Hybrid toolbox stuff?
rmpath(genpath('/home/akshatd/mine/umich/sem2/740/toolboxes/hybtbx-linux/'))

% plot MPT

figure();
sgtitle("MPT: x_1, x_2, x_3");

for i = 1:3
  subplot(3, 1, i);
  plot(times, data.X2(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.X20(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("x_" + num2str(i));
  legend('Location', 'best');
end

snapnow;

figure();
sgtitle("MPT: x_4, x_5, x_6");

for i = 4:6
  subplot(3, 1, i - 3);
  plot(times, data.X2(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.X20(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("x_" + num2str(i));
  legend('Location', 'best');
end

snapnow;

figure();
sgtitle("MPT: u_1, u_2, u_3");

for i = 1:3
  subplot(3, 1, i);
  plot(times, data.U2(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.U20(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("u_" + num2str(i));
  ylim([-0.12, 0.12]);
  legend('Location', 'best');
end

snapnow;

figure();
sgtitle('MPT: Computation Time');
plot(times, data.execTime2, 'DisplayName', 'N=2');
grid on;
hold on;
plot(times, data.execTime20, 'DisplayName', 'N=20');
xlabel('t [s]');
ylabel('Time [s]');
legend('Location', 'best');
snapnow;

% plot Hybrid

figure();
sgtitle("Hybrid Toolbox: x_1, x_2, x_3");

for i = 1:3
  subplot(3, 1, i);
  plot(times, data.X2h(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.X20h(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("x_" + num2str(i));
  legend('Location', 'best');
end

snapnow;

figure();
sgtitle("Hybrid Toolbox: x_4, x_5, x_6");

for i = 4:6
  subplot(3, 1, i - 3);
  plot(times, data.X2h(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.X20h(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("x_" + num2str(i));
  legend('Location', 'best');
end

snapnow;

figure();
sgtitle("Hybrid Toolbox: u_1, u_2, u_3");

for i = 1:3
  subplot(3, 1, i);
  plot(times, data.U2h(i, :), 'DisplayName', 'N=2');
  grid on;
  hold on;
  plot(times, data.U20h(i, :), 'DisplayName', 'N=20');
  xlabel('t [s]');
  ylabel("u_" + num2str(i));
  ylim([-0.12, 0.12]);
  legend('Location', 'best');
end

snapnow;

figure();
sgtitle('Hybrid Toolbox: Computation Time');
plot(times, data.execTime2h, 'DisplayName', 'N=2');
grid on;
hold on;
plot(times, data.execTime20h, 'DisplayName', 'N=20');
xlabel('t [s]');
ylabel('Time [s]');
legend('Location', 'best');
snapnow;

%% 2.a scdynamics function

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
