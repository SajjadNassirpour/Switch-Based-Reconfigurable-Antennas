clc
clear 
close all

% Number of users
K=4;
% Number of elements in the array antenna
M=64;
% SNR of the channel estimation process
SNR_est=20;
% SNR at the receiver
SNR_value=40;
% Numer of trials
iters=100;

% Optimization parameters
mu=10^(-4);
r=1e5;
divide_rate=2;

% noisy_cond=false : Noiseless channel are used in our algorithm  
% noisy_cond=true  : Noisy channel are used in our algorithm   
noisy_cond=false; 

% Generate the channels
[H_noisy,H_nonoise]=rayleigh_gen(K,M,iters,SNR_est);


for ii=1:length(SNR_value)
    ii
[sum_rate(:,ii)]=MMSE_app(K,M,iters,SNR_value(ii),H_noisy,H_nonoise)
[ref_final_rate(ii),ref_ratio]=Binary_ref_approach(K,M,iters,H_nonoise,SNR_value(ii),mu,r,divide_rate)
[our_final_rate(ii),our_ratio,TDMA_rate]=Binary_our_approach(K,M,iters,H_noisy,H_nonoise,SNR_value(ii),mu,r,divide_rate,noisy_cond)
end



