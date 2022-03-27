clc
clear 
close all

% Number of users
K=4;
% Number of elements in the array antenna
M=12;
% SNR of the channel estimation process
SNR_est=20;
% SNR at the receiver
SNR_value=-10:10:40;
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

% Optimization parameters of the genetic algorithm in the BIS method
Ns=2;
L=2;

% BF results with {0,1}
sum_rate_BF=zeros(1,length(SNR_value));

% BF results with {-1,0,1} in two timeslots
sum_rate_BF3_half=zeros(1,length(SNR_value));

% Our optimization method
our_final_rate=zeros(1,length(SNR_value));

% Optimization method with the filled function in [31]
ref1_final_rate=zeros(1,length(SNR_value));

% Optimization method with the filled function in [32]
ref2_final_rate=zeros(1,length(SNR_value));

% BIS optimization method with {0,1}
BIS_sum_rate=zeros(1,length(SNR_value));

% BIS optimization method with two timeslots and {-1,0,1}
BIS_sum_rate3=zeros(1,length(SNR_value));


for ii=1:length(SNR_value)
    
sum_rate_BF(ii)=BF_approach(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond);

sum_rate_BF3_half(ii)=BF3_half_approach(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond);

[ref1_final_rate(ii),ref1_complexity]=Binary_ref_approach(K,M,iters,H_nonoise,SNR_value(ii),mu,r,divide_rate);


[BIS_sum_rate3(ii), Switch_mat_BIS3]=BIS_low_complexity_3_half(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,Ns,L,noisy_cond);


[BIS_sum_rate(ii), Switch_mat_BIS]=BIS_low_complexity(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,Ns,L,noisy_cond);

[our_final_rate(ii),our_complexity,TDMA_rate]=Binary_our_approach(K,M,iters,H_noisy,H_nonoise,SNR_value(ii),mu,r,divide_rate,noisy_cond);

divide_rate=1.5;
[ref2_final_rate(ii),ref2_complexity]=Binary_ref2_approach(K,M,iters,H_nonoise,SNR_value(ii),mu,r,divide_rate);

end
sum_rate_BF
sum_rate_BF3_half
BIS_sum_rate3
BIS_sum_rate
our_final_rate
ref1_final_rate
ref2_final_rate