%% AE740 HW1 akshatdy

clc;
clear;
close all;

%% 1. Proof
%%
% $$J(x) = x^TQx + c^Tx$$
%
% Derivative of a scalar with respect to a vector
%
% $$\nabla_x c^Tx = c$$
%
% Product Rule
%
% $$\nabla_x x^TQx = (Q + Q^T)x$$
%
% since we know that Q is symmetric
%
% $$\nabla_x x^TQx = 2Qx$$
%
% therefore,
%
% $$\nabla_x J(x) = 2Qx + c$$

%% 1.a Second term
%%
% $$c^Tx = \sum_{i=1}^{n} c_i x_i$$
%
% partial derivatives with respect to $x_i$
%
% $$\frac{\partial}{\partial x_i}c^Tx = c_i$$
%
% therefore,
%
% $$\nabla_x c^Tx = c$$

%% 1.b First term
%%
% $$x^TQx = \sum_{i=1}^{n}\sum_{j=1}^{n} x_i Q_{ij} x_j$$
%
% partial derivatives with respect to $x_i$
%
% $$\frac{\partial}{\partial x_i}x^TQx = \sum_{j=1}^{n} Q_{ij} x_j + \sum_{j=1}^{n} Q_{ji} x_j$$
%
% since Q is symmetric, $Q_{ij} = Q_{ji}$
%
% $$\frac{\partial}{\partial x_i}x^TQx = \sum_{j=1}^{n} Q_{ij} x_j + \sum_{j=1}^{n} Q_{ij} x_j$$
%
% $$\frac{\partial}{\partial x_i}x^TQx = 2\sum_{j=1}^{n} Q_{ij} x_j$$
%
% therefore,
%
% $$\nabla_x x^TQx = 2Qx$$
%
% for example, when n = 2
%
% $$x^TQx = \left[\matrix{x_1 & x_2}\right]\left[\matrix{Q_{11} & Q_{12} \cr Q_{21} & Q_{22}}\right]\left[\matrix{x_1 \cr x_2}\right]$$
%
% $$x^TQx = x_1Q_{11}x_1 + x_1Q_{12}x_2 + x_2Q_{21}x_1 + x_2Q_{22}x_2$$
%
% since Q is symmetric, $Q_{12} = Q_{21}$
%
% $$x^TQx = x_1Q_{11}x_1 + 2x_1Q_{12}x_2 + x_2Q_{22}x_2$$
%
% partial derivatives with respect to $x_1$
%
% $$\frac{\partial}{\partial x_1}x^TQx = 2x_1Q_{11} + 2x_2Q_{12}$$
%
% partial derivatives with respect to $x_2$
%
% $$\frac{\partial}{\partial x_2}x^TQx = 2x_1Q_{12} + 2x_2Q_{22}$$
%
% therefore,
%
% $$\nabla_x x^TQx = \left[\matrix{2Q_{11}x_1 + 2Q_{12}x_2 \cr 2Q_{12}x_1 + 2Q_{22}x_2}\right]$$
%
% and
%
% $$2Qx = 2\left[\matrix{Q_{11} & Q_{12} \cr Q_{21} & Q_{22}}\right]\left[\matrix{x_1 \cr x_2}\right] = \left[\matrix{2Q_{11}x_1 + 2Q_{12}x_2 \cr 2Q_{21}x_1 + 2Q_{22}x_2}\right]$$
%
% therefore,
%
% $$\nabla_x x^TQx = 2Qx$$

%% 1.c Hessian
%%
% $$\nabla_x J(x) = 2Qx + c$$
%
% using the same rule as we used for the second part of the derivative
%
% $$\nabla_x^2 J(x) = 2Q$$
%