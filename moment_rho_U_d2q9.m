function [Rho,U] = moment_rho_U_d2q9(f)
%% Function [Rho, U] = moment_rho_U_d2q2(f) computes the physical density and velocity based on the given PDF at any location using D2Q9 lattice
%% f is the PDF and must be a row vector
%% Rho is the density and will be a scalar
%% U is the velocity and will be a column vector


% Define lattice velocity Xi (E)
E = [0, 1, 0, -1, 0, 1, -1, -1, 1;...
     0, 0, 1, 0, -1, 1, 1, -1, -1];

% Find Rho
Rho = sum(f);

% Find U
U = E*f'/Rho;