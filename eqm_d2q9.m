function f_eq = eqm_d2q9(Rho,U)
%% Function f_eq = eqm_d2q9(Rho,U) computes the equalibrium PDF based on the given Density and Velocity at any location using D2Q9 lattice
%% Rho is the density and must be a scalar
%% U is the velocity and must be a column vector
%% f_eq is the equalibrium PDF and will be a row vector

% Define latteral velocity xi (E), weight (w), speed of sound (c_s)
E = [0  0;
     1  0;
     0  1;
    -1  0;
     0 -1;
     1  1;
    -1  1;
    -1 -1;
     1 -1];
w = [4/9; 1/9; 1/9; 1/9; 1/9; 1/36; 1/36; 1/36; 1/36];
% c_s = 1/sqrt(3);

f_eq = (w.*Rho.*(1 + (3*E*U) + 9*(E*U).^2/2 - 3*(U(1)^2 + U(2)^2)/2))';