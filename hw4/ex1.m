clc;
clear;
close all;
clear variables
format shortG;

%% 1 Quadractic programming solver

%% 1.a Standard dual projected gradient algorithm

%%
% Equation in the question is
%
% $3x_1^2 + x_2^2 + 2x_1x_2 + x_1 + 6x_2 + 2$
%
% the 2 can be ignored from the optimization problem as it is a constant
% the resulting matrices that produce this equation are

H = [6 2;
     2 2];
q = [1; 6];

%%
% constraints are
%
% $2x_1 + 3x_2 \geq 4$
%
% $x_1 \geq 0$
%
% $x_2 \geq 0$
%
% have to reverse the sign of the inequality to make it  $\leq$

A = [-2 -3;
     -1 0;
     0 -1];
b = [-4; 0; 0];

% run with MATLAB's quadprog
[x, fval, exitflag, output, lam] = quadprog(H, q, A, b);
disp('1.a Standard dual projected gradient algorithm')
disp('MATLAB quadprog:');
disp('Solution =');
disp(x);
disp('Lagrange multipliers =');
disp(lam.ineqlin);

[x_my, lam_my] = myQP(H, q, A, b, zeros(size(A, 1), 1));
disp('myQP:');
disp('Solution =');
disp(x_my);
disp('Lagrange multipliers =');
disp(lam_my);

%% 1.c Discretize the mass-spring-damper system
disp('1.c Discrete time model of the mass-spring system')
% have to re-declare A and B with consants for m and k
Ts = 0.1;
m = 1;
k = 1;
Ac = [0 1;
      -k / m 0];
Bc = [0;
      1 / m];

[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, [], [], Ts, 'zoh');
disp('Ad =');
disp(Ad);
disp('Bd =');
disp(Bd);

%% 1.d Tracking MPC formulation
% $\Delta x_{1, k+1}$ and $\Delta x_{2, k+1}$ can be represented as $x_{k+1}$
%
% Then $\Delta x_{k+1} = x_{k+2} - x_{k+1}$
%
% $= Ax_{k+1} + Bu_{k+1} - (Ax_{k} + Bu_{k})$
%
% $= A(x_{k+1} - x_{k}) + B(u_{k+1} - u_{k})$
%
% Using definitions of $\Delta x_k$ and $\Delta u_k$
%
% $\Delta x_{k+1} = A\Delta x_k + B\Delta u_k$
%
% The error $e$ can be calculated by:
%
% $x_{k+1} = x_k + \Delta x_k$ -- (1)
%
% $e_k = x_k - r$ -- (2)
%
% $e_{k+1} = x_{k+1} - r$ -- (3)
%
% Substituting (1) into (3)
%
% $e_{k+1} = x_k + \Delta x_k - r$
%
% Using (2) to reduce this equation
%
% Thus: $e_{k+1} = e_k + \Delta x_k$
%
% $u_{k+1} = u_k + \Delta u_k$, so it just uses the $A_d$ to pick $u$ and $B_d$ to pick $\Delta u$.
%
% $x_{1, k+1} = x_{1, k} + \Delta x_{1, k}$ by just rearranging the given equation for $\Delta x_k$

A = [Ad, zeros(2, 3);
     1, 0, 1, 0, 0;
     0, 0, 0, 1, 0;
     1, 0, 0, 0, 1];
B = [Bd; 0; 1; 0];

%% 1.e LQ-MPC Problem formulation
Q = diag([0, 0, 1, 0, 0]);
R = 1;
[K, Pdxu, E] = dlqr(A(1:3, 1:3), B(1:3), Q(1:3, 1:3), R);
P = blkdiag(Pdxu, zeros(2, 2));
disp('1.e LQ-MPC Problem formulation')
disp('P =');
disp(P);

%% 1.f Represent constraints
lN = 10;
% large number
xlim.max = [lN, lN, lN, 0.2, 0.2]; % [dx1, dx2, e, u, x_1]
xlim.min = -xlim.max;
% Note our "control" is "control increment", actual control is the fourth state
ulim.max = [lN];
ulim.min = -ulim.max;

%% 1.h Simulate MPC closed loop

%% 1.a Standard dual projected gradient Function
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

  while k <= Nit % see Note 2
    lam = max(lam - 1 / L * df, 0);
    df = Hd * lam + qd;
    k = k + 1;
  end

  U = -invH * (G' * lam + q);
  return;
  % Note 1: Can pre-compute these quantities and pass as arguments to
  % myQP. In LQ-MPC setting, only q and Wtilde depend on x_0, makes no
  % sense to constantly re-compute these
  % Note 2: Change to another criterion as desired
end

%% 1.b Mass-spring dynamics
function xdot = msd(t, x, u, m_msd, k_msd)
  xdot = [
          x(2);
          (-k_msd / m_msd) * x(1) + (u / m_msd);
          ];
end

%% 1.g function that forms matrices needed for QP
function [H, L, G, W, T, IMPC] = formQPMatrices(A, B, Q, R, P, xlim, ulim, N)
  % This function forms the matrices needed for the constrained QP
  % Inputs:
  %   A, B: state-space matrices
  %   Q, R, P: cost function matrices
  %   xlim, ulim: state and input constraints
  %   N: prediction horizon

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

  H = S' * Qbar * S + Rbar;
  L = S' * Qbar * M;
  G = [S;
       -S;
       eye(N * nu);
       -eye(N * nu)];
  W = [
       repmat(xlim.max, N, 1);
       repmat(-xlim.min, N, 1);
       repmat(ulim.max, N, 1);
       repmat(-ulim.min, N, 1);
       ];
  T = [
       M;
       -M;
       zeros(nx, N * nu);
       zeros(nx, N * nu);
       ];
  IMPC = [eye(nu, nu), zeros(nu, (N - 1) * nu)];

end
