clc;
clear;
close all;

% setup matrices for part 2
A = [
	4/3 -2/3;
	1 0;
];
B = [
	1;
	0;
];
C = [-2/3 1];

% weighing matrices
Q = C' * C;
R = 0.001;
P = 0;

%% 2.a
disp("2.a");
if all(abs(eig(A)) < 1)
	disp('all eigenvalues of A are inside the unit circle');
else
	disp('all eigenvalues of A are not inside the unit circle');
end
%%
% A is a Schur matrix if all eigenvalues are inside the unit circle, so A is a Schur matrix

if rank(ctrb(A, B)) == size(A, 1)
	disp('(A, B) has full rank');
else
	disp('(A, B) does not have full rank');
end
%%
% (A, B) controllability matrix has full rank, so there are no uncontrollable modes, so it is stabilizable
%
% (A, B) is controllable if the controllability matrix has full rank, so it is controllable

if rank(obsv(A, C)) == size(A, 1)
	disp('(A, C) has full rank');
else
	disp('(A, C) does not have full rank');
end
%%
% (C, A) observability matrix has full rank, so there are no unobservable modes, so it is detectable
%
% (C, A) is observable if the observability matrix has full rank, so it is observable

%% 2.b Transfer function symbolically
disp("2.b");
z = sym('z');
disp("Symbolic transfer function");
disp(simplify(C * inv(z * eye(size(A)) - A) * B + 0));
disp("Eigenvalues of A");
disp(eig(A));
disp("Magnitude of eigenvalues of A");
disp(abs(eig(A)));
%%
% Transfer function is
%
% $$-\frac{2z - 3}{3z^2 - 4z + 2}$$
%
% Zeros(zeros of numerator)
% 
% $$2z - 3 = 0$$
%
% $$z = \frac{3}{2}$$
%
% $$z = 1.5$$
%
% Zeroes are out of the unit disk
%
% Poles(zeros of denominator)
%
% $$3z^2 - 4z + 2 = 0$$
%
% $$z = \frac{2}{3} \pm \frac{\sqrt{2}}{3}$$
%
% $$z = 0.6667 \pm 0.4714i$$
%
% Poles are inside the unit disk
%
% Poles are the same as eigenvalues of A

%% 2.c Step response
disp("2.c");
dstep(A, B, C, 0);
