% AE740 HW1 akshatdy
clc;
close all;

%% 2
Ac = [
	0.0000, 1.0000, 0.0000, 0.0000, 0.0000;
	0.0000, -0.86939, 43.2230, -17.2510, -1.5766;
	0.0000, 0.99335, -1.3411, -0.16897, -0.25183;
	0.0000, 0.0000, 0.0000, -20.0000, 0.0000;
	0.0000, 0.0000, 0.0000, 0.0000, -20.0000;
];
Bc = [
	0.0000, 0.0000;
	0.0000, 0.0000;
	0.0000, 0.0000;
	20.0000, 0.0000;
	0.0000, 20.0000;
];
Cc = [
	1.0000, 0.0000, 0.0000, 0.0000, 0.0000;
	1.0000, 0.0000, -1.0000, 0.0000, 0.0000;
];

%% 2.a check if model is open loop stable
eigAc = eig(Ac);
fprintf('2.a\nmodel is ');
% check if all eigenvalues are in left half plane
notLeftHalfPlane = eigAc(real(eigAc) >= 0);
if size(notLeftHalfPlane) ~= 0
	fprintf('open loop unstable, not all eigenvalues in left half plane:');
	display(notLeftHalfPlane);
else
	fprintf('open loop stable\n');
end

%% 2.b check if model is controllable
% check if rank of controllability matrix is equal to dimension of x
fprintf('2.b\ncontrollability matrix rank: %d\n', rank(ctrb(Ac, Bc)));
fprintf('dimension of x: %d\n', size(Ac, 1));
fprintf('model is ');
if rank(ctrb(Ac, Bc)) == size(Ac, 1)
	fprintf('controllable\n');
else
	fprintf('not controllable\n');
end
fprintf('time horizon is infinitely small since this is a continuous time system and input is unconstrained\n')

%% 2.c convert model to discrete time
Ts = 0.01;
Dc = zeros(size(Cc, 1), size(Bc, 2));
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, Ts, 'zoh');
fprintf('\n2.c');
display(Ad);
display(Bd);

%% 2.d convert to discrete time using expressions given in class
Ad2 = expm(Ac * Ts);
Bd2 = inv(Ac) * (Ad2 - eye(size(Ad2))) * Bc;
fprintf('\n2.d');
display(Ad2);
display(Bd2);
fprintf('Ac is a singular matrix, so it is not invertible. Hence we cannot use it to find Bd using the expressions in class\n');

%% 2.e check if discrete time model is open loop stable
eigAd = eig(Ad);
fprintf('\n2.e\ndiscrete time model is ');
% check if all eigenvalues are in unit circle
outsideUnitCircle = eigAd(abs(eigAd) >= 1);
if size(outsideUnitCircle) ~= 0
	fprintf('open loop unstable, eigenvalues outside unit circle:');
	display(outsideUnitCircle);
else
	fprintf('open loop stable\n');
end

%% 2.f check if discrete time model is controllable
fprintf('\n2.f\ncontrollability matrix rank: %d\n', rank(ctrb(Ad, Bd)));
fprintf('dimension of x: %d\n', size(Ad, 1));
fprintf('discrete time model is ');
if rank(ctrb(Ad, Bd)) == size(Ad, 1)
	fprintf('controllable\n');
else
	fprintf('not controllable\n');
end
fprintf('since n=%d, the time horizon is n*Ts=%.2fs\n', size(Ad, 1), size(Ad, 1) * Ts);

%% 2.g check if continuous time model is closed-loop stable
Fc = [
	-2.8900, 0.7780;
	1.9800, 3.3400;
];
Kc = [
	2.1100, 0.8906, 4.9107, -0.5343, -0.1009;
	-5.3200, -0.8980, -4.6618, 0.4280, 0.1099;
];

fprintf('\n2.g\n');
%%
% $\dot{x_c} = A_c x_c + B_c u_c$ where $u_c = F_c r + K_c x_c$
% 
% $$\dot{x_c} = A_c x_c + B_c (F_c r + K_c x_c)$$
% 
% $$\dot{x_c} = A_c x_c + B_c F_c r + B_c K_c x_c$$
% 
% $$\dot{x_c} = (A_c + B_c K_c) x_c + B_c F_c r$$
%
% check matrix $(A_c + B_c K_c)$
eigAcBcKc = eig(Ac + Bc * Kc);
notLeftHalfPlane = eigAcBcKc(real(eigAcBcKc) >= 0);
fprintf('system with feedforward and feedback is open loop ')
if size(notLeftHalfPlane) ~= 0
	fprintf('unstable, eigenvalues not in left half plane:');
	display(notLeftHalfPlane);
else
	fprintf('stable, all eigenvalues lie in the left half plane\n');
end

%%
% to get steady state gain, we set $\dot{x_c} = 0$
%
% $$0 = (A_c + B_c K_c) x_c + B_c F_c r$$
%
% $$x_c = -(A_c + B_c K_c)^{-1} B_c F_c r$$
%
% $$y_c = C_c x_c$$
%
% $$y_c = C_c (-(A_c + B_c K_c)^{-1} B_c F_c r)$$
%
% if $y_c = Hr$, then
%
% $$H = - C_c (A_c + B_c K_c)^{-1} B_c F_c$$

H = - Cc * inv(Ac + Bc * Kc) * Bc * Fc;
fprintf('steady state gain H:\n');
display(H);

%% 2.h continuous time closed-loop simulation for 5 seconds
fprintf('\n2.h\n');
tspan = [0 5];
x0 = zeros(size(Ac, 1), 1);

%% case 1
r = [0.01745; -0.01745];
odefun = @(t, x) f16modelCont(x, uSaturated(r, x, Fc, Kc), Ac, Bc);
[T, X] = ode45(odefun, tspan, x0);
Y = (Cc * X')';
U = uSaturated(r, X', Fc, Kc)';
figure(1);
hold on;
plot(T, Y);
plot(T, U);
title('2.h time history of states for case 1');
xlabel('time (s)');
ylabel('angle (rad)');
legend('pitch angle', 'flight path angle', 'elevator position', 'flaperon position');
hold off;

fprintf('Case 1: the closed loop system is stable, ');
fprintf('because the feedforward input is small enough ');
fprintf('to not saturate the total input to the system\n');

%%
% Case 1: the closed loop system is stable, 
% because the feedforward input is small enough 
% to not saturate the total input to the system
%%

%% case 2
r = [0.1745; -0.1745];
odefun = @(t, x) f16modelCont(x, uSaturated(r, x, Fc, Kc), Ac, Bc);
[T, X] = ode45(odefun, tspan, x0);
Y = (Cc * X')';
U = uSaturated(r, X', Fc, Kc)';
figure(2);
hold on;
plot(T, Y);
plot(T, U);
title('2.h time history of states for case 2');
xlabel('time (s)');
ylabel('angle (rad)');
ylim([-1 1]);
legend('pitch angle', 'flight path angle', 'elevator position', 'flaperon position');
hold off;

fprintf('Case 2: the closed loop system is not stable, ');
fprintf('because the feedforward inputs saturate the overall input, ');
fprintf('and the system essentially has no way to stabilize itself\n');

%% 2.i check if discrete time closed-loop model is stable
eigAdBdKc = eig(Ad + Bd * Kc);
outsideUnitCircle = eigAdBdKc(abs(eigAdBdKc) >= 1);
fprintf('\n2.i\ndiscrete time closed loop system is ');
if size(outsideUnitCircle) ~= 0
	fprintf('unstable, eigenvalues outside unit circle:');
	display(outsideUnitCircle);
else
	fprintf('stable, all eigenvalues lie in the unit circle\n');
end

%% 2.j check if discrete time closed-loop model is stable with Ts=0.5s
Ts5 = 0.5;
[Ad5, Bd5, Cd5, Dd5] = c2dm(Ac, Bc, Cc, Dc, Ts5, 'zoh');
eigAdBdKc = eig(Ad5 + Bd5 * Kc);
outsideUnitCircle = eigAdBdKc(abs(eigAdBdKc) >= 1);
fprintf('\n2.j with Ts=0.5, discrete time closed-loop model is ');
if size(outsideUnitCircle) ~= 0
	fprintf('unstable, eigenvalues outside unit circle:');
	display(outsideUnitCircle);
else
	fprintf('stable, all eigenvalues lie in the unit circle\n');
end

%% 2.k discrete time closed-loop simulation for 5 seconds
fprintf('\n2.k\n');
tspan = [0:Ts:5];

x0 = zeros(size(Ad, 1), 1);
% case 1
r = [0.01745; -0.01745];
X = zeros(length(tspan), length(x0));
X(1, :) = x0;
for t=2:length(tspan)
	X(t, :) = f16modelDisc(X(t-1, :)', r, Ad, Bd, Kc, Fc)';
end
Y = (Cd * X')';
U = uSaturated(r, X', Fc, Kc)';
figure(3);
hold on;
plot(tspan, Y);
plot(tspan, U);
title('2.k time history of states for discrete time model');
xlabel('time (s)');
ylabel('angle (rad)');
legend('pitch angle', 'flight path angle', 'elevator position', 'flaperon position');
hold off;

fprintf('The discrete time closed loop system is stable');

%% functions for all the parts, needs to be at the end

% 2.h
function xdot = f16modelCont(x, u, A, B)
	xdot = A*x + B*u;
end

function u = uSaturated(r, x, Fc, Kc)
	uMax = [0.4363; 0.3491];
	uMin = -uMax;
	u = Fc * r + Kc * x;
	u = max(uMin, min(uMax, u));
end
	
% 2.k
function xnext = f16modelDisc(x, r, Ad, Bd, Kc, Fc)
	xnext = Ad*x + Bd*uSaturated(r, x, Fc, Kc);
end

%% ?