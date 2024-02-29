clc;
clear;
close all;
clear variables
format shortG;

%% 3 Control of a car

%% 3.a Discrete time model
disp('3.a Discrete time model');
% consts
m = 1891;
vx = 20;
Cf = -18e3;
Cr = -28e3;
lf = 1.5;
lr = 1.55;
Iz = 3200;
nx = 2;
nu = 1;

Ac = [2 * (Cf + Cr) / (m * vx) 2 * (Cf * lf - Cr * lr) / (m * vx) - vx;
      2 * (lf * Cf - lr * Cr) / (Iz * vx) 2 * (Cf * lf^2 + Cr * lr^2) / (Iz * vx)];
Bc = [-2 * Cf / m;
      -2 * lf * Cf / Iz];
Cc = [1 / vx, lf / vx;
      1 / vx, -lr / vx;
      ];
Dc = [
      -1;
      0;
      ];
Ts = 20e-3;
card = c2d(ss(Ac, Bc, Cc, Dc), Ts, 'zoh');
Ad = card.A;
disp('Ad = ');
disp(Ad);
Bd = card.B;
disp('Bd = ');
disp(Bd);
Cd = card.C;
Dd = card.D;

% check if Ad is schur
if all(abs(eig(Ad)) < 1)
  disp('Ad is schur');
else
  disp('Ad is not schur');
end

if rank(ctrb(Ad, Bd)) == size(Ad, 1)
  disp('(Ad, Bd) has full rank, it is controllable');
else
  disp('(Ad, Bd) does not have full rank, it is not controllable');
end

%% 3.b simulate using MPT
disp('3.b simulate using MPT');
Ax = [Ad, Bd, 0 * Ad(:, 1);
      0 * Ad(1, :), 1, 0;
      0 * Ad(1, :), 0, 1;
      ];
Bx = [Bd;
      1;
      0;
      ];
% Cx is af and ar with stuff for u and r
Cx = [Cd, Dd, 0 * Cd(:, 1); ];
Dx = zeros(size(Cx, 1), size(Bx, 2));
Ex = [0, 1, 0, -1];
Qy = 3;
R = 0.1;
Q = Ex' * Qy * Ex;
P = zeros(size(Ax, 1));

N = 10; % Horizon

% MPT stuff, do once
tbxmanager restorepath
mpt_init
% create the model
modelMPT = LTISystem('A', Ax, 'B', Bx, 'C', Cx, 'D', Dx, 'Ts', Ts);
modelMPT.x.min = [-3, -1, -0.15, -inf];
modelMPT.x.max = -modelMPT.x.min;
modelMPT.u.min = -0.06;
modelMPT.u.max = -modelMPT.u.min;
modelMPT.y.min = [-0.1, -0.07];
modelMPT.y.max = -modelMPT.y.min;
modelMPT.u.penalty = QuadFunction(R);
modelMPT.x.penalty = QuadFunction(Q);
modelMPT.x.with('terminalPenalty');
modelMPT.x.terminalPenalty = QuadFunction(P);
ctrl = MPCController(modelMPT, N);

%% 3.c Simulate using Hybrid toolbox
disp('3.c Simulate using Hybrid toolbox');
%%
% Both the toolboxes give the same results, so we can use either one for
% simulation, and it assures us that the simulation is correct.
% We also see that the constraints are satisfied.

% Hybrid toolbox stuff
addpath( ...
  genpath('/home/akshatd/mine/umich/sem2/740/toolboxes/hybtbx-linux/'), ...
'-end')
clear cost constraints horizon;
% have to make a new C matrix to pass back x
Cxh = [Cx;
       eye(4);
       ];
Dxh = zeros(size(Cxh, 1), size(Bx, 2));

cost.Q = Q;
cost.R = R;
cost.P = P;
cost.rho = Inf;
% constraints.xmin = [-3, -1, -0.15, -inf];
% constraints.xmax = -constraints.xmin;
constraints.umin = -0.06;
constraints.umax = -constraints.umin;
% add the constraints for x here
constraints.ymin = [-0.1, -0.07, -3, -1, -0.15, -inf];
constraints.ymax = -constraints.ymin;

horizon.Nu = N;
horizon.N = N;
horizon.Ncu = N - 1;
horizon.Ncy = N - 1;
cardx = ss(Ax, Bx, Cxh, Dxh, Ts);
ctrlh = lincon(cardx, 'reg', cost, horizon, constraints, 'qpact', 0);

x0 = [0; 0; 0; 0.2];
y = Cx * x0;
yh = Cxh * x0;
Tsim = 2;
times = 0:Ts:Tsim;

data.X = zeros(4, Tsim / Ts + 1);
data.X(:, 1) = x0;
data.Xh = zeros(4, Tsim / Ts + 1);
data.Xh(:, 1) = x0;
data.Y = zeros(2, Tsim / Ts + 1);
data.Y(:, 1) = y;
data.Yh = zeros(6, Tsim / Ts + 1);
data.Yh(:, 1) = yh;

x = x0;
xh = x0;
dataIdx = 2;

for t = times(1: end - 1)

if t < 1
  x(4) = 0.2;
  xh(4) = 0.2;
else
  x(4) = -0.2;
  xh(4) = -0.2;
end

[u, feasible, openloop] = ctrl.evaluate(x);
% data.U(:, dataIdx) = u;
x = Ax * x + Bx * u;
data.X(:, dataIdx) = x;
data.Y(:, dataIdx) = Cx * x;
% calc side slip angle
data.X(1, dataIdx) = x(1) / vx;

uh = eval(ctrlh, xh);
xh = Ax * xh + Bx * uh;
data.Xh(:, dataIdx) = xh;
data.Yh(:, dataIdx) = Cxh * xh;
% calc side slip angle
data.Xh(1, dataIdx) = xh(1) / vx;

dataIdx = dataIdx + 1;
end

% some other Hybrid toolbox stuff?
rmpath(genpath('/home/akshatd/mine/umich/sem2/740/toolboxes/hybtbx-linux/'))

% plot MPT
figure();
sgtitle("MPT states");

subplot(3, 1, 1);
plot(times, data.X(2, :), "DisplayName", "Yaw rate");
grid on;
hold on;
plot(times, data.X(4, :), "DisplayName", "Reference yaw rate");
xlabel('t [s]');
ylabel("Yaw rate");
ylim([-0.3, 0.3]);
legend("location", "best");

subplot(3, 1, 2);
plot(times, data.X(1, :));
grid on;
xlabel('t [s]');
ylabel("Side slip angle");

subplot(3, 1, 3);
plot(times, data.X(3, :));
grid on;
xlabel('t [s]');
ylabel("Steering angle");

snapnow;

figure();
sgtitle("MPT outputs");

subplot(2, 1, 1);
plot(times, data.Y(1, :));
grid on;
xlabel('t [s]');
ylabel("Front slip angle");
ylim([-0.12, 0.12]);

subplot(2, 1, 2);
plot(times, data.Y(2, :));
grid on;
xlabel('t [s]');
ylabel("Rear slip angle");
ylim([-0.07, 0.07]);

snapnow;

% plot Hybrid toolbox

figure();
sgtitle("Hybrid MPC states");

subplot(3, 1, 1);
plot(times, data.Xh(2, :), "DisplayName", "Yaw rate");
grid on;
hold on;
plot(times, data.Xh(4, :), "DisplayName", "Reference yaw rate");
xlabel('t [s]');
ylabel("Yaw rate");
ylim([-0.3, 0.3]);
legend("location", "best");

subplot(3, 1, 2);
plot(times, data.Xh(1, :));
grid on;
xlabel('t [s]');
ylabel("Side slip angle");

subplot(3, 1, 3);
plot(times, data.Xh(3, :));
grid on;
xlabel('t [s]');
ylabel("Steering angle");

snapnow;

figure();
sgtitle("Hybrid MPC outputs");

subplot(2, 1, 1);
plot(times, data.Yh(1, :));
grid on;
xlabel('t [s]');
ylabel("Front slip angle");
ylim([-0.12, 0.12]);

subplot(2, 1, 2);
plot(times, data.Yh(2, :));
grid on;
xlabel('t [s]');
ylabel("Rear slip angle");
ylim([-0.07, 0.07]);

snapnow;
