clc
clear 
close all

K=4;  % Number of users (BTSs)
M=64; % NUmber of elements in array antenna

est_SNR_mat=20;  % SNR(dB) in channel estimation process

iters=100;  % Number of iterations

pow_mat=-10:10:40;  % SNR interval

% Optimization parameters in our approach
mu=0.0001;
r=10^5;
divide_rate=2;

% Sum-rate with Analog beamforming
sum_rate_analog=zeros(length(est_SNR_mat),1,length(pow_mat));  

% Sum-rate with TDMA+Analog beamforming
sum_rate_tdma_analog=zeros(length(est_SNR_mat),1,length(pow_mat));

% Sum-rate with Beamnulling+Analog beamforming
sum_rate_beam_null=zeros(length(est_SNR_mat),1,length(pow_mat));

% Sum-rate with MMSE
sum_rate_mmse=zeros(length(est_SNR_mat),2,length(pow_mat));

% Sum-rate with our approach
sum_rate_our=zeros(length(est_SNR_mat),1,length(pow_mat));



for jj=1:length(est_SNR_mat)
    est_SNR=est_SNR_mat(jj);
    
    % Generate the channels
    [H_noisy,H_nonoise,beam_vector, null_vector]=real_world_channel_beam_null(K,M,iters,est_SNR, 'uniform');

    % MMSE
    sum_rate_mmse(jj,:,:)=MMSE_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise);
    mmse_iter=squeeze(sum_rate_mmse(jj,:,:))
    
    % Analog Beamforming
    sum_rate_analog(jj,:,:)=Analog_app_complex(K,M,iters,pow_mat,H_nonoise,beam_vector);
    analog_iter=squeeze(sum_rate_analog(jj,:,:))'
    
    % TDMA + Analog Beamforming
    sum_rate_tdma_analog(jj,:,:)=TDMA_Analog_app_complex(K,M,iters,pow_mat,H_nonoise,beam_vector);
    tdma_analog_iter=squeeze(sum_rate_tdma_analog(jj,:,:))'
    
    % Beamnulling + Analog Beamforming
    sum_rate_beam_null(jj,:,:)=Beam_null_app_complex(K,M,iters,pow_mat,H_nonoise,null_vector);
    beam_null_iter=squeeze(sum_rate_beam_null(jj,:,:))'

    % noise_cond
    % false: Noiseless channels
    % true: Noisy channels 
    noise_cond=true; 
    
    % Switch_mat : Optimal switch configuration 
    [sum_rate_our(jj,:,:), Switch_mat_R]=Binary_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise,mu,r,divide_rate,noise_cond);
    our_iter=squeeze(sum_rate_our(jj,:,:))'
end

