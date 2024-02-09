clc;
clear;
close all;

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

%% 3.b Compute feedback gain using uncMPC
disp("3.b")
NHoriz = 20;
specRadList = zeros(20, 1);
eigvalList = zeros(NHoriz*2, 1);
K20 = []; % for use in 3.d
KList = zeros(NHoriz, 2); % for use in 3.e
for N = 1:NHoriz
	[S, M, Qbar, Rbar, K0N] = uncMPC(N, A, B, Q, R, P);
	KList(N, :) = K0N;
	eigval = eig(A + B * K0N);
	eigvalRow = (N-1) * 2 + 1;
	eigvalList(eigvalRow:eigvalRow+1, :) = eigval;
	specRad = max(abs(eigval));
	specRadList(N, 1) = specRad;
	% check if the max of abs of eigenvalues of A + BK are less than 1
	if specRad < 1
		disp("System is closed-loop stable at N = " + num2str(N));
	end
	if N == NHoriz
		K20 = K0N;
	end
end
%%
% The closed loop system is stable at N = 1 and N >= 11

% plot the spectral radius vs N
figure(1);
plot(1:NHoriz, specRadList);
grid on;
title('3.b Spectral Radius vs N');
xlabel('N');
ylabel('Spectral Radius');
ylim([0.6, 1.6]);

%% 3.c Plot eigenvalues in complex plane for N[1, 20]
figure(2);
hold on;
for N = 1:NHoriz
	eigvalRow = (N-1) * 2 + 1;
	plot(real(eigvalList(eigvalRow:eigvalRow+1, :)), imag(eigvalList(eigvalRow:eigvalRow+1, :)), "*");
end
grid on;
title('3.c Eigenvalues in Complex Plane');
xlabel('Real');
xlim([-0.3, 1.7]);
ylabel('Imaginary');
ylim([-0.6, 0.6]);
% plot unit circle
th = 0:0.01:2*pi;
x = cos(th);
y = sin(th);
plot(x, y, 'r');
legend("N = 1", "N = 2", "N = 3", "N = 4", "N = 5", "N = 6", "N = 7", "N = 8", "N = 9", "N = 10", "N = 11", "N = 12", "N = 13", "N = 14", "N = 15", "N = 16", "N = 17", "N = 18", "N = 19", "N = 20", "Unit Circle", "Location", "eastoutside");
hold off;

%% 3.d Generate LQR gain using dlqr
disp("3.d")
[Kinf, P, e] = dlqr(A, B, Q, R, P);
display("Kinf = " + mat2str(-Kinf));
display("K20 = " + mat2str(K20));
%%
% Yes, $K_{0, N}$ is close to $K_{\infty}$ and $K_{0, N} \rightarrow K_{\infty}$ as $N \rightarrow \infty$ which means it is stabilizing.

%% 3.e Plot the components of K
figure(3);
hold on;
grid on;
title('3.e K(1) vs N');
xlabel('N');
ylabel('K(1)');
plot(2:NHoriz, KList(2:end, 1));
% plot the limits as the Kinf values
plot([2, NHoriz], [-Kinf(1), -Kinf(1)], '--r');
legend("KN", "Kinf", "Location", "best")
hold off;

figure(4);
hold on;
grid on;
title('3.e K(2) vs N');
xlabel('N');
ylabel('K(2)');
plot(2:NHoriz, KList(2:end, 2));
% plot the limits as the Kinf values
plot([2, NHoriz], [-Kinf(2), -Kinf(2)], '--r');
legend("KN", "Kinf", "Location", "best")

%%
% Plotting without N = 1 to display the graph better
%
% Yes, $K_{0, N}$ converges to $K_{\infty}$ as $N \rightarrow \infty$

%% 3.a Implement unconstrained LQ-MPC
function [S, M, Qbar, Rbar, K0N] = uncMPC(N, A, B, Q, R, P)
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

		% Compute the first column of S
		for i = 1:N
			rowStart = (i - 1) * nx + 1;
			rowEnd = i * nx;
			S(rowStart:rowEnd, 1:nu) = A^(i-1)*B;
		end
		% Pad the first column and set it to other columns of S
		for i = 2:N
			% rowStart = (i - 1) * nx + 1;
			% rowEnd = i * nx;
			colStart = (i - 1) * nu + 1;
			colEnd = i * nu;
			zeroRows = (i - 1) * nx;
			zeroCols = nu;
			S(:, colStart:colEnd) = [zeros(zeroRows, zeroCols); S(1:end - zeroRows, 1:nu)];
		end

		% Compute first row of M
		M(1:nx, :) = A;
		% Compute the rest of M
		for i = 2:N
			rowStart = (i - 1) * nx + 1;
			rowEnd = i * nx;
			% just multiply the previous rows by A to get higher powers
			M(rowStart:rowEnd, :) = A * M(rowStart - nx:rowEnd - nx, :);
		end

		% Compute Qbar except for the last row
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

		% Compute Rbar
		for i = 1:N
			% R is square so we can reuse indices
			rowStart = (i - 1) * nu + 1;
			rowEnd = i * nu;
			Rbar(rowStart:rowEnd, rowStart:rowEnd) = R;
		end

		% Compute K0N
		K0N = -[eye(nu), zeros(nu, nu * (N - 1))] * inv(S'*Qbar*S + Rbar) * S'*Qbar*M;
end
