% AE740 HW1 akshatdy
clc;
close all;

%% 3
n = 0.00114;
Ac = [
	0 0	0	1	0 0;
	0 0	0	0	1 0;
	0 0 0 0 0 1;
	3*n^2 0	0	0	2*n 0;
	0 0	0	-2*n 0 0;
	0 0	-n^2 0	0 0;
];
Cc = [
	0.5 0.5 0 0 0 0;
	0 0 1 0 0 0;
];

%% 3.a get Ad, Cd
Ts = 30;
[Ad, Bd, Cd, Dd] = c2dm(Ac, [], Cc, [], Ts, 'zoh');
fprintf('3.a\n');
display(Ad);
display(Cd);

%% 3.b check if model is observable
fprintf('3.b\nobservability matrix rank: %d\n', rank(obsv(Ad, Cd)));
fprintf('number of states: %d\n', size(Ad, 1));
if rank(obsv(Ad, Cd)) == size(Ad, 1)
	fprintf('system is observable\n');
else
	fprintf('system is not observable\n');
end

%% 3.c design Luenberger observer
E = [
	0.9712;
	0.9377;
	0.9348 + 0.0686i;
	0.9348 - 0.0686i;
	0.9471 + 0.0753i;
	0.9471 - 0.0753i;
];
Ld = -place(Ad', Cd', E)';
fprintf('3.c\n');
display(Ld);

%% 3.d simulate discrete time
Xd = [0.1; 1; 2; 0;	0.01;	-0.01];
Xdhat = [0; 0; 0; 0; 0; 0];

Traj = [];
for k = 0:1:100
	Traj.time(k+1) = k;
	Traj.Xd(k+1,:) = Xd;
	Traj.Xdhat(k+1,:) = Xdhat;
	Yd = Cd*Xd;	% true output measurement
	Ydhat = Cd*Xdhat; % estimated output measurement
	Xd = Ad*Xd; % update model
	Xdhat = Ad*Xdhat + Ld*(Ydhat - Yd); % update the observer
end

figure(1);
hold on;
plot(Traj.time, Traj.Xd(:,2));
plot(Traj.time, Traj.Xdhat(:,2));
title('3.d: discrete time system for the nominal and observer (Y)');
xlabel('time (s)');
ylabel('Y (km)');
legend('true Y', 'estimated Y');
hold off;

figure(2);
hold on;
plot(Traj.time, Traj.Xd(:,5));
plot(Traj.time, Traj.Xdhat(:,5));
title('3.d: discrete time system for the nominal and observer (Ydot)');
xlabel('time (s)');
ylabel('Ydot (km/s)');
legend('true Ydot', 'estimated Ydot');
hold off;

figure(3);
hold on;
plot(Traj.time, Traj.Xd(:,3));
plot(Traj.time, Traj.Xdhat(:,3));
title('3.d: discrete time system for the nominal and observer (Z)');
xlabel('time (s)');
ylabel('Z (km)');
legend('true Z', 'estimated Z');
hold off;

figure(4);
hold on;
plot(Traj.time, Traj.Xd(:,6));
plot(Traj.time, Traj.Xdhat(:,6));
title('3.d: discrete time system for the nominal and observer (Zdot)');
xlabel('time (s)');
ylabel('Zdot (km/s)');
legend('true Zdot', 'estimated Zdot');
hold off;

%% 3.e simulate discrete time with noise
Xd = [0.1; 1; 2; 0;	0.01;	-0.01];
Xdhat = [0; 0; 0; 0; 0; 0];

Traj = [];
for k = 0:1:100
	Traj.time(k+1) = k;
	Traj.Xd(k+1,:) = Xd;
	Traj.Xdhat(k+1,:) = Xdhat;
	Yd = Cd*Xd+ (randn(2,1)*0.5);	% true output measurement with noise
	Ydhat = Cd*Xdhat; % estimated output measurement
	Xd = Ad*Xd; % update model
	Xdhat = Ad*Xdhat + Ld*(Ydhat - Yd); % update the observer
end

figure(5);
hold on;
plot(Traj.time, Traj.Xd(:,2));
plot(Traj.time, Traj.Xdhat(:,2));
title('3.e: noisy system for the nominal and observer (Y)');
xlabel('time (s)');
ylabel('Y (km)');
legend('true Y', 'estimated Y');
hold off;

figure(6);
hold on;
plot(Traj.time, Traj.Xd(:,5));
plot(Traj.time, Traj.Xdhat(:,5));
title('3.e: noisy system for the nominal and observer (Ydot)');
xlabel('time (s)');
ylabel('Ydot (km/s)');
legend('true Ydot', 'estimated Ydot');
hold off;

figure(7);
hold on;
plot(Traj.time, Traj.Xd(:,3));
plot(Traj.time, Traj.Xdhat(:,3));
title('3.e: noisy system for the nominal and observer (Z)');
xlabel('time (s)');
ylabel('Z (km)');
legend('true Z', 'estimated Z');
hold off;

figure(8);
hold on;
plot(Traj.time, Traj.Xd(:,6));
plot(Traj.time, Traj.Xdhat(:,6));
title('3.e: noisy system for the nominal and observer (Zdot)');
xlabel('time (s)');
ylabel('Zdot (km/s)');
legend('true Zdot', 'estimated Zdot');
hold off;

