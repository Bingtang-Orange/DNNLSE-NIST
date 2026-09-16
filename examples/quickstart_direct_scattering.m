% QUICKSTART_DIRECT_SCATTERING
% Solves the direct scattering problem for a single spectral parameter lambda
% and computes the scattering matrix and reflection coefficients.

clear; clc; close all;
addpath(genpath('../src'));
warning('off', 'chebfun:compose');

%% 1. Configuration and Parameters
L      = 40.0;       % Domain truncation half-length
lambda = 0.5;        % Spectral parameter (real or complex)

% Benchmark potential q(x) (Gaussian profile)
q0 = @(x) exp(-x.^2);

fprintf('=====================================================\n');
fprintf(' Direct Scattering Quickstart Demo at lambda = %.4f\n', lambda);
fprintf('=====================================================\n');

%% 2. Compute Reflection Coefficients via Chebfun
[rho, rho_hat] = compute_reflection(lambda, q0, L);

%% 3. Display Results
fprintf('Input potential : q0(x) = exp(-x^2)\n');
fprintf('Domain length   : [-%.1f, %.1f]\n', L, L);
fprintf('-----------------------------------------------------\n');
fprintf('rho     = %+.8f %+.8fi  (|rho|     = %.8f)\n', real(rho), imag(rho), abs(rho));
fprintf('rho_hat = %+.8f %+.8fi  (|rho_hat| = %.8f)\n', real(rho_hat), imag(rho_hat), abs(rho_hat));
fprintf('=====================================================\n');