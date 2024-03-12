clc;
clear;
close all;
clear variables
format shortG;

%% 2 Soft Constraints

%% 2.b Simulate MPC with soft constraints

% Define the base system
Ts = 0.1;
m = 1;
k = 1;
Ac = [0 1;
      -k / m 0];
Bc = [0;
      1 / m];

% Discretize the system
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, [], [], Ts, 'zoh');

% construct the extended system for the LQ-MPC problem
A = [Ad, zeros(2, 3);
     1, 0, 1, 0, 0;
     0, 0, 0, 1, 0;
     1, 0, 0, 0, 1];
B = [Bd; 0; 1; 0];

Q = diag([0, 0, 1, 0, 0]);
R = 1;
[K, Pdxu, E] = dlqr(A(1:3, 1:3), B(1:3), Q(1:3, 1:3), R);
P = blkdiag(Pdxu, zeros(2, 2));

cmds = [-0.19 0.19 -0.25];
cmd = cmds(1);
prev_cmd = cmd;
xk_msd_sm = [0; 0];
xk_msd_lg = [0; 0];
xk_sm = [0; 0; -cmd; 0; 0];
xk_lg = [0; 0; -cmd; 0; 0];

N = 15;
lN = 10;
% large number
XLIM.max = [lN, lN, lN, 0.2, 0.2]; % [dx1, dx2, e, u, x_1]
XLIM.min = -XLIM.max;
% Note our "control" is "control increment", actual control is the fourth state
ulim.max = [lN];
ulim.min = -ulim.max;
[H_sm, L_sm, G_sm, W_sm, T_sm, IMPC_sm] = formQPMatrices(A, B, Q, R, P, XLIM, ulim, N, 1e-3);
[H_lg, L_lg, G_lg, W_lg, T_lg, IMPC_lg] = formQPMatrices(A, B, Q, R, P, XLIM, ulim, N, 1e3);

lam_sm = ones(size(G_sm, 1), 1);
lam_lg = ones(size(G_lg, 1), 1);

Tcmd = 12;
Tsim = Tcmd * length(cmds);
times = 0:Ts:Tsim;
fidelity = Ts / 10;

data.x_msd_sm = zeros(2, Tsim / Ts);
data.x_msd_lg = zeros(2, Tsim / Ts);
data.u_sm = zeros(1, Tsim / Ts);
data.u_lg = zeros(1, Tsim / Ts);
data.r = zeros(1, Tsim / Ts);
data_idx = 1;

for t = times
  % get the current command
  cmd = cmds(min(1 + floor(t / Tcmd), length(cmds)));

  if cmd ~= prev_cmd
    % update the error in the state
    xk_sm(3) = xk_sm(3) - cmd + prev_cmd;
    xk_lg(3) = xk_lg(3) - cmd + prev_cmd;
  end

  % solve the QP
  [U_sm, lam_sm] = myQP(H_sm, L_sm * xk_sm, G_sm, W_sm + T_sm * xk_sm, lam_sm);
  [U_lg, lam_lg] = myQP(H_lg, L_lg * xk_lg, G_lg, W_lg + T_lg * xk_lg, lam_lg);

  % get the first control increment
  delta_uk_sm = IMPC_sm * U_sm;
  delta_uk_lg = IMPC_lg * U_lg;
  uk_sm = xk_sm(4) + delta_uk_sm; % nu is just 1 so this works
  uk_lg = xk_lg(4) + delta_uk_lg; % nu is just 1 so this works

  % simulate the system
  [~, xk1_msd_ode_sm] = ode45(@(t, x) msd(t, x, uk_sm, 1, 1), [t + fidelity:fidelity:t + Ts], xk_msd_sm);
  [~, xk1_msd_ode_lg] = ode45(@(t, x) msd(t, x, uk_lg, 1, 1), [t + fidelity:fidelity:t + Ts], xk_msd_lg);
  xk1_msd_sm = xk1_msd_ode_sm(end, :)';
  xk1_msd_lg = xk1_msd_ode_lg(end, :)';

  % update the state
  delta_xk_msd_sm = xk1_msd_sm - xk_msd_sm;
  delta_xk_msd_lg = xk1_msd_lg - xk_msd_lg;

  xk1_sm = [delta_xk_msd_sm; % delta x1 and x2
            xk_sm(3) + delta_xk_msd_sm(1); % e
            uk_sm; % u
            xk1_msd_sm(1)]; % x1

  xk1_lg = [delta_xk_msd_lg; % delta x1 and x2
            xk_lg(3) + delta_xk_msd_lg(1); % e
            uk_lg; % u
            xk1_msd_lg(1)]; % x1

  % update
  xk_msd_sm = xk1_msd_sm;
  xk_msd_lg = xk1_msd_lg;
  xk_sm = xk1_sm;
  xk_lg = xk1_lg;
  prev_cmd = cmd;

  % save data
  data.x_msd_sm(:, data_idx) = xk_msd_sm;
  data.x_msd_lg(:, data_idx) = xk_msd_lg;
  data.u_sm(data_idx) = uk_sm;
  data.u_lg(data_idx) = uk_lg;
  data.r(data_idx) = cmd;

  data_idx = data_idx + 1;

end

% plot the results
fig = figure();
fig.Position(3:4) = [800, 600];
sgtitle('2.b Simulate MPC with slack variable');
subplot(3, 1, 1);
plot(times, data.x_msd_sm(1, :), 'DisplayName', '$\mu = 1e-3$', "LineWidth", 2);
grid on;
hold on;
plot(times, data.x_msd_lg(1, :), 'Color', '#77AC30', 'DisplayName', '$\mu = 1e3$', "LineWidth", 2);
plot(times, data.r, '--k', 'DisplayName', 'reference', "LineWidth", 1);
yline(0.2, '--r', 'DisplayName', '$x_1$ max', "LineWidth", 1);
yline(-0.2, '--r', 'DisplayName', '$x_1$ min', "LineWidth", 1);
xlabel('t [s]');
ylabel('x_1');
xlim tight;
ylim([-0.3, 0.3]);
legend("Location", "best", "Interpreter", "latex");

subplot(3, 1, 2);
plot(times, data.x_msd_sm(2, :), 'DisplayName', '$\mu = 1e-3$', "LineWidth", 2);
grid on;
hold on;
plot(times, data.x_msd_lg(2, :), 'Color', '#77AC30', 'DisplayName', '$\mu = 1e3$', "LineWidth", 2);
xlabel('t [s]');
ylabel('x_2');
xlim tight;
ylim([-0.5, 0.5]);
legend("Location", "best", "Interpreter", "latex");

subplot(3, 1, 3);
stairs(times, data.u_sm, 'DisplayName', '$\mu = 1e-3$', "LineWidth", 2);
grid on;
hold on;
stairs(times, data.u_lg, 'Color', '#77AC30', 'DisplayName', '$\mu = 1e3$', "LineWidth", 2);
xlabel('t [s]');
ylabel('u');
yline(0.2, '--r', 'DisplayName', 'u max', "LineWidth", 1);
yline(-0.2, '--r', 'DisplayName', 'u min', "LineWidth", 1);
xlim tight;
ylim([-0.6, 0.6]);
legend("Location", "best", "Interpreter", "latex");

disp('1.h Simulate MPC closed loop');
fprintf('x_1_sm limits = [%f  %f]\n', min(data.x_msd_sm(1, :)), max(data.x_msd_sm(1, :)));
fprintf('x_1_lg limits = [%f  %f]\n', min(data.x_msd_lg(1, :)), max(data.x_msd_lg(1, :)));
fprintf('u_sm limits = [%f  %f]\n', min(data.u_sm), max(data.u_sm));
fprintf('u_lg limits = [%f  %f]\n', min(data.u_lg), max(data.u_lg));

%%
% With the slack variables added, we can see that the constraints for both values of
% $\mu$ are violated. This makes sense since the slack variables are added to the cost function
% to allow for some constraint violation. We can also see that the constraints are violated more
% for the smaller value of $\mu$ than for the larger value of $\mu$. This is expected since the
% larger value of $\mu$ penalizes the slack variables more, and thus the cost function will be
% more sensitive to the constraint violations.
%
% We can also observe that $x_1$ reaches its reference value faster for the smaller value of $\mu$
% than for the larger value of $\mu$ because the controller is able to be more aggressive with the
% smaller value of $\mu$ by violating the constraints more on the controller input(u). We can also
% observe that the final referrence value which violates the limits for $x_1$ is actually reached
% when using the smaller value of $\mu$. The controller is able to keep violating the constraints
% on the controller input(u) and $x_1% as it balances the cost of following the reference with
% the cost of violating the constraints on the controller input(u) and $x_1%.

%% 2.a function that forms matrices needed for QP with soft constraints
function [H, L, Gs, W, T, IMPC] = ...
    formQPMatrices(A, B, Q, R, P, xlim, ulim, N, slackPenalty)
  % This function forms the matrices needed for the constrained QP
  % Inputs:
  %   A, B: state-space matrices
  %   Q, R, P: cost function matrices
  %   xlim, ulim: state and input constraints
  %   N: prediction horizon
  %   slackPenalty: penalty for slack variables

  nx = size(A, 1);
  nu = size(B, 2);

  S = zeros(N * nx, N * nu);
  % Compute the first column of S
  for i = 1:N
    rowStart = (i - 1) * nx + 1;
    rowEnd = i * nx;
    S(rowStart:rowEnd, 1:nu) = A^(i - 1) * B;
  end

  % Pad the first column and set it to other columns of S
  for i = 2:N
    colStart = (i - 1) * nu + 1;
    colEnd = i * nu;
    zeroRows = (i - 1) * nx;
    zeroCols = nu;
    S(:, colStart:colEnd) = [zeros(zeroRows, zeroCols); S(1:end - zeroRows, 1:nu)];
  end

  % add padding to account for the slack variables
  Ss = [S, zeros(N * nx, nx)];

  M = zeros(N * nx, nx);
  % Compute first row of M
  M(1:nx, :) = A;
  % Compute the rest of M
  for i = 2:N
    rowStart = (i - 1) * nx + 1;
    rowEnd = i * nx;
    % just multiply the previous rows by A to get higher powers
    M(rowStart:rowEnd, :) = A * M(rowStart - nx:rowEnd - nx, :);
  end

  Qbar = zeros(N * nx, N * nx);
  % Compute Qbar and set the last row to P
  for i = 1:N
    % Q is square so we can reuse indices
    rowStart = (i - 1) * nx + 1;
    rowEnd = i * nx;
    temp = Q;

    if i == N
      temp = P;
    end

    Qbar(rowStart:rowEnd, rowStart:rowEnd) = temp;
  end

  Rbar = zeros(N * nu, N * nu);
  % Compute Rbar
  for i = 1:N
    % R is square so we can reuse indices
    rowStart = (i - 1) * nu + 1;
    rowEnd = i * nu;
    Rbar(rowStart:rowEnd, rowStart:rowEnd) = R;
  end

  % add penalty for slack variables
  Rbars = [Rbar, zeros(N * nu, nx);
           zeros(nx, N * nu), slackPenalty * eye(nx)];

  H = Ss' * Qbar * Ss + Rbars;
  L = Ss' * Qbar * M;
  Gs = [S, repmat(eye(nx), N, 1); % eye should be nx * ne, but theyre the same
        -S, repmat(eye(nx), N, 1);
        eye(N * nu), zeros(N * nu, nx);
        -eye(N * nu), zeros(N * nu, nx);
        ];

  W = [
       repmat(xlim.max', N, 1);
       repmat(-xlim.min', N, 1);
       repmat(ulim.max', N, 1);
       repmat(-ulim.min', N, 1);
       ];
  T = [
       -M;
       M;
       zeros(N * nu, nx);
       zeros(N * nu, nx);
       ];
  IMPC = [eye(nu, nu), zeros(nu, (N - 1) * nu), zeros(nu, nx)];

end

% Standard dual projected gradient Function
function [U, lam] = myQP(H, q, A, b, lam0)
  % This function implements the dual projected
  % gradient algorithm for solving a QP problem.
  % Minimize 1/2 * U' * H * U + q' * U subject to G * U <= Wtilde

  % compared to his notes, G=A, U=x, W=b
  G = A; Wtilde = b;
  invH = inv(H); G_invH = G * invH; % see Note 1
  Hd = G_invH * G';
  % see Note 1
  qd = G_invH * q + Wtilde;
  Nit = 30;
  % maximum number of iterations
  lam = lam0;
  L = norm(Hd);
  k = 1;
  df = Hd * lam + qd;

  while k <= Nit
    lam = max(lam - 1 / L * df, 0);
    df = Hd * lam + qd;
    k = k + 1;
  end

  U = -invH * (G' * lam + q);
  return;
end

% Mass-spring dynamics
function xdot = msd(t, x, u, m_msd, k_msd)
  xdot = [
          x(2);
          (-k_msd / m_msd) * x(1) + (u / m_msd);
          ];
end
