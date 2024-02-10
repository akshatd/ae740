clc;
clear;
close all;

%% 4

% setup matrices from part 2
A = [
	4/3 -2/3;
	1 0;
];
B = [
	1;
	0;
];
C = [-2/3 1];
Q = C' * C;
R = 0.001;
P = 0;

%% 4.b Compare using problem 2
disp("4.b")
[Kinf, Pinf, e] = dlqr(A, B, Q, R, P);
disp("Using dlqr: Pinf=");disp(Pinf);
% for n=[2, 11, 20]
% 	[S, M, Qbar, Rbar, K0N, P0N] = uncMPCric(n, A, B, Q, R, P);
% 	disp("Using uncMPCric: N= " + num2str(n) + " P= " + mat2str(P0N) + " norm(Pinf-P0N)= " + num2str(norm(Pinf-P0N)));
% end
% cant do loops cos then matlab wont publish the disp
disp("Using uncMPCric:");
n=2;
[S, M, Qbar, Rbar, K02, P0N] = uncMPCric(n, A, B, Q, R, P);
disp("N=" + num2str(n));disp("P=");disp(P0N);disp("norm(Pinf-P0N)=");disp(norm(Pinf-P0N));
n=11;
[S, M, Qbar, Rbar, K011, P0N] = uncMPCric(n, A, B, Q, R, P);
disp("N=" + num2str(n));disp("P=");disp(P0N);disp("norm(Pinf-P0N)=");disp(norm(Pinf-P0N));
n=20;
[S, M, Qbar, Rbar, K020, P0N] = uncMPCric(n, A, B, Q, R, P);
disp("N=" + num2str(n));disp("P=");disp(P0N);disp("norm(Pinf-P0N)=");disp(norm(Pinf-P0N));
%%
% The difference in the norm between P from dlqr and uncMPCric decreases as N increases and is 0.00020802 at N=20, which means that the two methods are converging.

%% 4.c Compare using P = Pinf
disp("4.c")
disp("Using dlqr: Kinf=");disp(Kinf);disp("Pinf=");disp(Pinf);
% for n=[2, 11, 20]
% 	[S, M, Qbar, Rbar, K0N, P0N] = uncMPCric(n, A, B, Q, R, Pinf);
% 	disp("Using uncMPCric: N= " + num2str(n) + " K0N= " + mat2str(K0N) + " norm(Kinf-K0N)= " + num2str(norm(-Kinf-K0N)) + " P= " + mat2str(P0N) + " norm(Pinf-P0N)= " + num2str(norm(Pinf-P0N)));
% end
% cant do loops cos then matlab wont publish the disp
disp("Using uncMPCric:");
n=2;
[S, M, Qbar, Rbar, K0N, P0N] = uncMPCric(n, A, B, Q, R, Pinf);
disp("N=" + num2str(n));disp("K0N=");disp(K0N);disp("norm(Kinf-K0N)= ");disp(norm(-Kinf-K0N));disp("P=");disp(P0N);disp("norm(Pinf-P0N)=");disp(norm(Pinf-P0N));
n=11;
[S, M, Qbar, Rbar, K0N, P0N] = uncMPCric(n, A, B, Q, R, Pinf);
disp("N=" + num2str(n));disp("K0N=");disp(K0N);disp("norm(Kinf-K0N)= ");disp(norm(-Kinf-K0N));disp("P=");disp(P0N);disp("norm(Pinf-P0N)=");disp(norm(Pinf-P0N));
n=20;
[S, M, Qbar, Rbar, K0N, P0N] = uncMPCric(n, A, B, Q, R, Pinf);
disp("N=" + num2str(n));disp("K0N=");disp(K0N);disp("norm(Kinf-K0N)= ");disp(norm(-Kinf-K0N));disp("P=");disp(P0N);disp("norm(Pinf-P0N)=");disp(norm(Pinf-P0N));
%%
% The difference in the norm between K from dlqr and uncMPCric is extremely small regardless of N, in the order of 1e-15.
%
% The same is true is for P as well, so we can conclude that using P = Pinf in uncMPCric makes both P and K converge to the values from dlqr regardless of N.
%
% By setting P=Pinf, the gains K for any prediction horizon N are are stabilizing.

%% 4.d Proof
% $$J = \sum_{k=0}^{\infty} x_k^T Q x_k + u_k^T R u_k$$
%
% (Convert $u$ to $x$) We have $u_k = K_k x_k$, substituting into $J$:
%
% $$J = \sum_{k=0}^{\infty} x_k^T Q x_k + K^T x_k^T R K x_k$$
%
% $$J = \sum_{k=0}^{\infty} x_k^T (Q + K^T R K) x_k$$
%
% (Convert $x_k$ to $x_0$) For closed-loop system, given control gain $K$:
% 
% $$x_{k+1} = Ax_k + Bu_k$$
%
% $$x_{k+1} = Ax_k + BKx_k$$
%
% $$x_{k+1} = (A+BK)x_k$$
%
% $$x_{k+1} = (A+BK)x_k \rightarrow x_{k+1} = (A+BK)^{k+1} x_0 \rightarrow x_k = (A+BK)^k x_0$$
%
% substituting into $J$:
%
% $$J = \sum_{k=0}^{\infty} x_0^T ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k x_0$$
%
% $$J = x_0^T \left(\sum_{k=0}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k\right) x_0$$
%
% Let $P_k = \sum_{k=0}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$, then:
%
% Then $J = x_0^T P_k x_0$
%
% $$P_k = \sum_{k=0}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% Multiply by $(A+BK)^T$ on the left and $(A+BK)$ on the right:
%
% $$(A + BK)^T P_k (A + BK) = (A + BK)^T \left(\sum_{k=0}^{\infty} (A+BK)^k)^T (Q + K^T R K) (A+BK)^k\right) (A + BK)$$
%
% $$(A + BK)^T P_k (A + BK) = \sum_{k=0}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1}$$
%
% Subtract $P_k$ from both sides:
%
% $$(A + BK)^T P_k (A + BK) - P_k = \sum_{k=0}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1} - P_k$$
%
% $$= \sum_{k=0}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1} - \sum_{k=0}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% For $k=0$, $(A + BK)^T P_k (A + BK) - P_k$:
%
% $$= ((A+BK)^1)^T (Q + K^T R K) (A+BK)^1 - ((A+BK)^0)^T (Q + K^T R K) (A+BK)^0$$
%
% $$= (A+BK)^T (Q + K^T R K) (A+BK) - (Q + K^T R K)$$
%
% For $k=[1, \infty]$, $(A + BK)^T P_k (A + BK) - P_k$:
% 
% $$= \sum_{k=1}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1} - \sum_{k=1}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% Adding the two gives us:  $(A + BK)^T P_k (A + BK) - P_k$:
%
% $$= (A+BK)^T (Q + K^T R K) (A+BK) - (Q + K^T R K)$$
% 
% $$+\sum_{k=1}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1} - \sum_{k=1}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% $$= - (Q + K^T R K) $$
% 
% $$+(A+BK)^T (Q + K^T R K) (A+BK) + \sum_{k=1}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1}$$
% 
% $$- \sum_{k=1}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% We know that if $k=0$, the first term in the summation:
% 
% $$\sum_{k=0}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1} = (A+BK)^T (Q + K^T R K) (A+BK)$$
%
% So:
%
% $$(A+BK)^T (Q + K^T R K) (A+BK) + \sum_{k=1}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1}$$
%
% $$= \sum_{k=0}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1}$$
% 
% additionally,
%
% $$\sum_{k=0}^{\infty} ((A+BK)^{k+1})^T (Q + K^T R K) (A+BK)^{k+1} = \sum_{k=1}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% so, $(A + BK)^T P_k (A + BK) - P_k$:
%
% $$= - (Q + K^T R K) $$
% 
% $$+ \sum_{k=1}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
% 
% $$- \sum_{k=1}^{\infty} ((A+BK)^k)^T (Q + K^T R K) (A+BK)^k$$
%
% $$= - (Q + K^T R K)$$
%
% Hence:
%
% $$(A + BK)^T P_k (A + BK) - P_k = - (Q + K^T R K)$$ 
%
% the Lyapunov equation is:
%
% $$(A+BK)^T P_k (A+BK) - P_k + (Q + K^T R K) = 0$$
%
% Rewriting the equation:
%
% $$(A+BK)^T P_k (A+BK) - P_k = - (Q + K^T R K)$$
%
% Which is the same result we got earlier, hence, $P_k$ is the solution to the Lyapunov equation.

%% 4.e Compute finite horizon performance
disp("4.e")
x0 = [1;-0.5];
% For N=inf
Kinf = -Kinf;
ABK = A + B*Kinf;
QKRK = Q + Kinf'*R*Kinf;
Pk = dlyap(ABK', QKRK);
Jinf = x0'*Pk*x0;
disp("Cost under N=inf: " + num2str(Jinf));
% For N=11
ABK = A + B*K011;
QKRK = Q + K011'*R*K011;
Pk = dlyap(ABK', QKRK);
J11 = x0'*Pk*x0;
disp("Cost under N=11: " + num2str(J11));
%%
% The cost under N=inf is a lot lower than N=11, which is expected since the cost is minimized over a much longer horizon.

%% 4.a Riccati recursion
function [S, M, Qbar, Rbar, K0N, P0N] = uncMPCric(N, A, B, Q, R, P)
	% Compute the matrices S, M, Qbar, Rbar, and K0N
	% for the unconstrained LQ-MPC problem
	% 
	% Inputs:
	%   N: Prediction horizon
	%   A: State matrix
	%   B: Input matrix
	%   Q: State cost matrix
	%   R: Input cost matrix
	%   P: Terminal state cost matrix

	nx = size(A, 1);
	nu = size(B, 2);

	% Initialize matrices
	S = zeros(N*nx, N*nu);
	M = zeros(N*nx, nx);
	Qbar = zeros(N*nx, N*nx);
	Rbar = zeros(N*nu, N*nu);

	Pk = P;

	for k=1:N
		Kkm1 = -inv(R + B'*Pk*B)*B'*Pk*A;
		Pkm1 = Q + A'*Pk*A + A'*Pk*B*Kkm1;

		Pk = Pkm1;
	end
	K0N = Kkm1;
	P0N = Pkm1;
end
