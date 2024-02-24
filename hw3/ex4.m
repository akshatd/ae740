clc;
clear;
close all;
clear variables

%% 4 Spacecraft Motion

% consts
n = 0.00114;

%% 4.a Open loop dynamics
disp('4.a Open loop dynamics');
Ac = [0, 0, 0, 1, 0, 0;
      0, 0, 0, 0, 1 0;
      0, 0, 0, 0, 0, 1;
      3 * n^2, 0, 0, 0, 2 * n, 0;
      0, 0, 0, -2 * n, 0, 0;
      0, 0, -n^2, 0, 0, 0;
      ];
disp('Ac = ');
disp(Ac);

eigAc = eig(Ac);
% check if all eigenvalues are in left half plane
notLeftHalfPlane = eigAc(real(eigAc) >= 0);

if size(notLeftHalfPlane) ~= 0
  fprintf('model is open loop unstable, not all eigenvalues in left half plane:');
  display(notLeftHalfPlane);
else
  fprintf('model is open loop stable\n');
end

%% 4.b Observability
disp('4.b Observability')
Cc = [1, 0, 0, 0, 0, 0;
      0, 1, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0;
      ];
disp('C = ');
disp(Cc);

if rank(obsv(Ac, Cc)) == size(Ac, 1)
  disp('(C, A) has full rank, it is observable');
else
  disp('(C, A) does not have full rank, it is not observable');
end

%% 4.c Discrete time model
disp('4.c Discrete time model');

Ts = 30;
% Bc = zeros(size(Ac, 1), 1);
% Dc = zeros(size(Cc, 1), 1);

[Ad, ~, Cd, ~] = c2dm(Ac, [], Cc, [], Ts, 'zoh');
disp('Ad = ');
disp(Ad);
disp('Cd = ');
disp(Cd);

%% 4.d Moving Horizon Observer
disp('4.d Moving Horizon Observer');
P = diag([1, 1, 1, 0.001, 0.001, 0.001]);
R = diag([0.01, 0.01, 0.01]);
Q = diag([0.001, 0.001, 0.001, 0.1, 0.1, 0.1]);

X0 = [1; 1; -1; 0.002; -0.002; 0.004];
X0_hat = [0; 0; 0; 0; 0; 0; ];

X = X0;
X_hat = X0_hat;
Y = Cd * X + 0.01 * randn(3, 1);

data.X = zeros(6, 101);
data.X(:, 1) = X;
data.X_hat = zeros(6, 101);
data.X_hat(:, 1) = X_hat;

for k = 1:100
  % observe
  Y = Cd * X + 0.01 * randn(3, 1);

  % estimate
  X_hat_bar = Ad * X_hat;
  P_bar = Q + Ad * P * Ad';
  L = P_bar * Cd' * inv(Cd * P_bar * Cd' + R);
  X_hat = X_hat_bar + L * (Y - Cd * X_hat_bar);
  P = P_bar - P_bar * Cd' * inv(Cd * P_bar * Cd' + R) * Cd * P_bar;

  % propagate
  X = Ad * X;
  data.X(:, k + 1) = X;
  data.X_hat(:, k + 1) = X_hat;
end

% plot
figure;
plot(0:100, data.X(1, :), "DisplayName", "x");
grid on;
hold on;
plot(0:100, data.X(2, :), "DisplayName", "y");
plot(0:100, data.X(3, :), "DisplayName", "z");
plot(0:100, data.X_hat(1, :), "DisplayName", "$\hat{x}$");
plot(0:100, data.X_hat(2, :), "DisplayName", "$\hat{y}$");
plot(0:100, data.X_hat(3, :), "DisplayName", "$\hat{z}$");
xlabel('k');
ylabel('Relative Position [km]');
legend("Location", "best", "Interpreter", "latex", "FontSize", 14);
title('Position states');

figure;
plot(0:100, data.X(4, :), "DisplayName", "$\dot{x}$");
grid on;
hold on;
plot(0:100, data.X(5, :), "DisplayName", "$\dot{y}$");
plot(0:100, data.X(6, :), "DisplayName", "$\dot{z}$");
plot(0:100, data.X_hat(4, :), "DisplayName", "$\hat{\dot{x}}$");
plot(0:100, data.X_hat(5, :), "DisplayName", "$\hat{\dot{y}}$");
plot(0:100, data.X_hat(6, :), "DisplayName", "$\hat{\dot{z}}$");
xlabel('k');
ylabel('Relative Velocity [km]');
legend("Location", "best", "Interpreter", "latex", "FontSize", 14);
title('Velocity states');

% "$\hat{\dot{x}}$"
