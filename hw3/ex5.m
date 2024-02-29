%% 5 Spacecraft Dynamics

clc;
clear;
close all;
clear variables
format shortG;

% consts
n = 1.107e-3;
m = 100;

%% 5.a State space form
disp("5.a");
Ac = [
      0, 0, 1, 0;
      0, 0, 0, 1;
      3 * n^2, 0, 0, 2 * n;
      0, 0, -2 * n, 0;
      ];
Bc = [
      0, 0;
      0, 0;
      1 / m, 0;
      0, 1 / m;
      ];

disp("Ac = ");
disp(Ac);
disp("Bc = ");
disp(Bc)
eigAc = eig(Ac);
% check if all eigenvalues are in left half plane
notLeftHalfPlane = eigAc(real(eigAc) >= 0);

if size(notLeftHalfPlane) ~= 0
  fprintf('model is open loop unstable, not all eigenvalues in left half plane:');
  display(notLeftHalfPlane);
else
  fprintf('model is open loop stable\n');
end

if rank(ctrb(Ac, Bc)) == size(Ac, 1)
  disp('(Ac, Bc) has full rank, it is controllable');
else
  disp('(Ac, Bc) has full rank, it is not controllable');
end

%% 5.b Discrete time
disp("5.b");
Ts = 20;
[Ad, Bd, ~, ~] = c2dm(Ac, Bc, [], [], Ts, 'zoh');
disp("Ad = ");
disp(Ad);
disp("Bd = ");
disp(Bd);

%% 5.e Minimize fuel consumption cost using inputs
disp("5.e");
N = 15;
q = 10^4;
X0 = [0; 10; 0; 0];
nx = size(Ad, 1);
nu = size(Bd, 2);

% using var names from the question

% XN = P X0 + R U
% P = A^N
P = Ad^N;
% R = [A^(N-1) B, A^(N-2) B, ..., B]
% R = nx x (N)(nu)
R = zeros(nx, N * nu);

for i = 1:N
  R(:, (i - 1) * nu + 1:i * nu) = Ad^(N - i) * Bd;
end

% now Jtilde = (1/2)Z' GAMMA1 Z + GAMMA2  Z
% GAMMA1 = 2 q Rprime' Rprime
% = 2 q [R' R 0; 0 0]
GAMMA1 = 2 * q * ...
  [
 R' * R, zeros(N * nu, N * nu);
 zeros(N * nu, N * nu), zeros(N * nu, N * nu)
 ];

% GAMMA2 = 2 q X0' P' R ones(1, N*nu)
GAMMA2 = [2 * q * X0' * P' * R, ones(1, N * nu)];

% from the matlab docs, we minimize (1/2) x' H x + f' x in x = quadprog(H,f)
% where x = Z, H = GAMMA1, f = GAMMA2

% subject to Ax <= b, so
% cai >= U -> U <= cai -> U - cai <= 0
% cai >= -U -> -U <= cai -> -U - cai <= 0
% U's are in the first half of Z, cai's are in the second half of Z
% we need to set up the matrix A so that it minuses them off
A = [eye(N * nu), -eye(N * nu);
     -eye(N * nu), -eye(N * nu)];
b = zeros(2 * N * nu, 1);

Zopt = quadprog(GAMMA1, GAMMA2, A, b);
% get Uopt from Zopt in the format of nu * N
Uopt = reshape(Zopt(1:N * nu), nu, N);

data.X = zeros(nx, N + 1);
data.X(:, 1) = X0;

for k = 2:N
  data.X(:, k) = Ad * data.X(:, k - 1) + Bd * Uopt(:, k - 1);
end

% plot trajectory
figure();
plot(data.X(2, :), data.X(1, :));
grid on;
xlabel('y');
ylabel('x');
title('5.e Position trajectory of the spacecraft');

figure();
plot(data.X(4, :), data.X(3, :));
grid on;
xlabel('y dot');
ylabel('x dot');
title('5.e Velocity trajectory of the spacecraft');

% plot inputs
figure();
stairs(Ts * (1:N), Uopt(1, :), "DisplayName", "u_x");
grid on;
hold on;
stairs(Ts * (1:N), Uopt(2, :), "DisplayName", "u_y");
xlabel("Time (s)");
ylabel("Input");
xlim([0, Ts * N + 10]);
title("5.e Input trajectory of the spacecraft");
legend("Location", "best");

%% 5.f Fuel Cost
disp("5.f");
J = sum(abs(Uopt), "all");
disp("5.f Fuel cost = ");
disp(J);

%% 5.g Additional Constraints
disp("5.g");
% now we do quadprog(H,f,A,b,Aeq,beq,lb,ub), ignore Aeq, beq, lb
% we need to set up ub so that infnorm(U) <= ub
% which is just |U| <= ub and leave cai to be inf
ub = [0.05 * ones(N * nu, 1);
      inf * ones(N * nu, 1)];
lb = [-0.05 * ones(N * nu, 1);
      -inf * ones(N * nu, 1)];
Zopt = quadprog(GAMMA1, GAMMA2, A, b, [], [], lb, ub);
Uopt = reshape(Zopt(1:N * nu), nu, N);

data.X = zeros(nx, N + 1);
data.X(:, 1) = X0;

for k = 2:N
  data.X(:, k) = Ad * data.X(:, k - 1) + Bd * Uopt(:, k - 1);
end

% plot trajectory
figure();
plot(data.X(2, :), data.X(1, :));
grid on;
xlabel('y');
ylabel('x');
title('5.g Position trajectory of the spacecraft');

figure();
plot(data.X(4, :), data.X(3, :));
grid on;
xlabel('y dot');
ylabel('x dot');
title('5.g Velocity trajectory of the spacecraft');

% plot inputs
figure();
stairs(Ts * (1:N), Uopt(1, :), "DisplayName", "u_x");
grid on;
hold on;
stairs(Ts * (1:N), Uopt(2, :), "DisplayName", "u_y");
xlabel("Time (s)");
ylabel("Input");
xlim([0, Ts * N + 10]);
title("5.g Input trajectory of the spacecraft");
legend("Location", "best");

% Fuel Cost
J = sum(abs(Uopt), "all");
disp("5.g Fuel cost = ");
disp(J);

%%
% The fuel cost is higher with the additional constraints, as expected. The
% additional constraints limit the inputs, so the spacecraft can't move as
% freely. The spacecraft has to use additional fuel over more time to achieve
% the same target state, thus the higher fuel cost. The position and velocity
% trajectories are also different, as expected, as the spacecraft can't move
% as freely.
