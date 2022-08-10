clc
clear 
close all

% Number of users
K=4;
% Number of elements in the array antenna
M=6;
% SNR of the channel estimation process
SNR_est=20;
% SNR at the receiver
SNR_value=-10:10:40;
% Numer of trials
iters=100;

% Optimization parameters
% We use mu to denote \epsilon
mu=10^(-3);
r=1e0;
divide_rate=100;

% noisy_cond=false : Noiseless channel are used in our algorithm  
% noisy_cond=true  : Noisy channel are used in our algorithm   
noisy_cond=false; 

% Generate the channels
[H_noisy,H_nonoise]=rayleigh_gen(K,M,iters,SNR_est);

% Optimization parameters of the BIS method
Ns=2;
L=2;

% Brute force with {0,1}
sum_rate_BF=zeros(1,length(SNR_value));

% Brute force with {-1,0,1} in two timeslots
sum_rate_BF3_half=zeros(1,length(SNR_value));

% Our optimization method
our_final_rate=zeros(1,length(SNR_value));

% Filled function in [31]
ref1_final_rate=zeros(1,length(SNR_value));

% Filled function in [32]
ref2_final_rate=zeros(1,length(SNR_value));

% BIS with {0,1}
BIS_sum_rate=zeros(1,length(SNR_value));

% BIS with two timeslots and {-1,0,1}
BIS_sum_rate3=zeros(1,length(SNR_value));

% Genetic algorithm (GA)
% -- Optimization parameter -- %
pop_size=100;
num_point_cross=1;
mut_prob=0.1;
sum_rate_GA=zeros(1,length(SNR_value));

% Simplified exhaustive search (SES)
% -- Optimization parameter -- %
bunch_size=4;
sum_rate_simplified_BF=zeros(1,length(SNR_value));

% Successive refinemnet(SR)
% -- Optimization parameter -- %
threshold=10^(-3);
rand_ini_num=100;
sum_rate_SR=zeros(1,length(SNR_value));

for ii=1:length(SNR_value)
   
    % --- Brute force with {0,1} --- %
    sum_rate_BF(ii)=BF_approach(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond);
    
    % --- Brute force with {-1,0,1} in two timeslots --- %
    sum_rate_BF3_half(ii)=BF3_half_approach(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond);
    
    % --- GA --- %
    [sum_rate_GA(ii), GA_complexity]=Genetic_approach(K,M,pop_size,num_point_cross,mut_prob,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond);

    % --- SES --- %
    [sum_rate_simplified_BF(ii), sim_BF_complexity]=Simplified_BF_approach(K,M,bunch_size,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond);

    % --- SR --- %
    [sum_rate_SR(ii), complexity_SR]=SR_approach(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,noisy_cond,threshold,rand_ini_num);

    % --- BIS with two timeslots and {-1,0,1} --- %
    [BIS_sum_rate3(ii), Switch_mat_BIS3]=BIS_low_complexity_3_half(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,Ns,L,noisy_cond);

    % --- BIS with {0,1} --- %
    [BIS_sum_rate(ii), Switch_mat_BIS,BIS_complexity]=BIS_low_complexity(K,M,iters,SNR_value(ii),H_noisy,H_nonoise,Ns,L,noisy_cond);
    
    % --- Our optimization method --- %
    initial_rand_num=1;
    [our_final_rate(ii),our_complexity,TDMA_rate]=Binary_our_approach(K,M,iters,H_noisy,H_nonoise,SNR_value(ii),mu,r,divide_rate,noisy_cond,initial_rand_num);

    % --- Filled function in [31] --- %
    divide_rate=500;
    [ref1_final_rate(ii),ref1_complexity]=Binary_ref_approach(K,M,iters,H_nonoise,SNR_value(ii),mu,r,divide_rate);

    % --- Filled function in [32] --- %
    divide_rate=10;
    [ref2_final_rate(ii),ref2_complexity]=Binary_ref2_approach(K,M,iters,H_nonoise,SNR_value(ii),mu,r,divide_rate);
end
