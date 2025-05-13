function [spike_t, poissonP_t] = simVS(f,vs,mu,lambda_mean,srate,dur_s,nsweep)
% Models a Poisson phase locked spike train after Kessler et al. 2021.
% N.B. Outputs spike times in milliseconds for compatibility with data.

% f - modulation frequency
% t - time
t = [0:1/srate:dur_s];
nt = length(t);

% Concentration parameter.
k_max = 100; dk = 1e-1;

% Calculate Io_k
dx = 1e-3; x = -pi:dx:pi;
Io_k = @(k) (1/(2*pi))*sum(exp(  k * cos(x) ))*dx;  

% Calculate I1_k
I1_k = @(k) (1/(2*pi))*sum(exp(  k .*cos(x) ).*cos(x))*dx;  

% Compute the VS for corresponding values of k. 
k_table = [0:dk:k_max];
VS_table =  arrayfun(@(k) I1_k(k)/Io_k(k),k_table);  

% This is a close match for Figure 2 of Kessler et al. 2021. 
% figure; plot(k_table,VS_table); pause;

% Find the correct value of k to use to give desired vector strength.
k = k_table(findnearest(vs,VS_table));

% calculate von Mises density function.
Pkf_t =  exp(k * cos(2*pi*t*f - mu))*(f/Io_k(k));

% Calculate the spike rate of a correspoding Poisson process (spikes per second).
lambda_t = lambda_mean * Pkf_t /f;

% Probability function
poissonP_t = lambda_t/srate;

for s = 1:nsweep
    % Multiply by 1000 to put in milliseconds.  
    spike_t{s} =  1000*t( find( rand(1,nt) < poissonP_t ) ) ;
end;