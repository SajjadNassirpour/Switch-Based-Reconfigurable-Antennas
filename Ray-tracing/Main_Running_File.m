clc
clear 
close all

K=4;  % Number of users (BTSs)
M=64; % NUmber of elements in array antenna

est_SNR_mat=20;  % SNR(dB) in channel estimation process

iters=100;  % Number of iterations

pow_mat=-10:10:40;  % SNR interval

% ee=0 means ZF receiver
% ee=1 means MMSE receiver
ee=1;   


% Optimization parameters in our approach
mu=0.0001;
r=10^5;
divide_rate=2;

% Sum-rate with MMSE
sum_rate_mmse=zeros(length(est_SNR_mat),2,length(pow_mat));

% Sum-rate with lattice-beamforming
sum_rate_lattice_inf_bits=zeros(length(est_SNR_mat),2,length(pow_mat));
sum_rate_lattice_1_bit=zeros(length(est_SNR_mat),2,length(pow_mat));
sum_rate_lattice_2_bits=zeros(length(est_SNR_mat),2,length(pow_mat));

% Sum-rate with our WLLS
sum_rate_WLLS_inf_bits=zeros(length(est_SNR_mat),2,length(pow_mat));
sum_rate_WLLS_1_bit=zeros(length(est_SNR_mat),2,length(pow_mat));
sum_rate_WLLS_2_bits=zeros(length(est_SNR_mat),2,length(pow_mat));

% Sum-rate with our approach
sum_rate_our=zeros(length(est_SNR_mat),2,length(pow_mat));


for jj=1:length(est_SNR_mat)
    est_SNR=est_SNR_mat(jj);
    
    % Generate the channels
    [H_noisy,H_nonoise,beam_vector, null_vector]=real_world_channel_beam_null(K,M,iters,est_SNR, 'uniform');

    % MMSE
    sum_rate_mmse(jj,:,:)=MMSE_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise);
    mmse_iter=squeeze(sum_rate_mmse(jj,:,:))
    
    
    
    % --- Inf - bits --- %
    bit=0;  % Number of bits in phase shifters (infinite bits)
    
    %-- Lattice --$
    sum_rate_lattice_inf_bits(jj,:,:)=HBF_Lattice(K,M,iters,pow_mat,H_noisy,H_nonoise,bit,ee);
    sum_rate_iter_lattice_inf_bits=squeeze(sum_rate_lattice_inf_bits(jj,:,:))
    
    %-- WLLS --$
    sum_rate_WLLS_inf_bits(jj,:,:)=HBF_WLLS(K,M,iters,pow_mat,H_noisy,H_nonoise,bit);
    sum_rate_iter_WLLS_inf_bits=squeeze(sum_rate_WLLS_inf_bits(jj,:,:))
    
    
    % --- 1 - bit --- %
    bit=1;
    
    %-- Lattice --$
    sum_rate_lattice_1_bit(jj,:,:)=HBF_Lattice(K,M,iters,pow_mat,H_noisy,H_nonoise,bit,ee);
    sum_rate_iter_lattice_1_bit=squeeze(sum_rate_lattice_1_bit(jj,:,:))
    
    %-- WLLS --$
    sum_rate_WLLS_1_bit(jj,:,:)=HBF_WLLS(K,M,iters,pow_mat,H_noisy,H_nonoise,bit);
    sum_rate_iter_WLLS_1_bit=squeeze(sum_rate_WLLS_1_bit(jj,:,:))

    
    
    % --- 2 - bits --- %
    bit=2;
    
    %-- Lattice --$
    sum_rate_lattice_2_bits(jj,:,:)=HBF_Lattice(K,M,iters,pow_mat,H_noisy,H_nonoise,bit,ee);
    sum_rate_iter_lattice_2_bits=squeeze(sum_rate_lattice_2_bits(jj,:,:))
    
    %-- WLLS --$
    sum_rate_WLLS_2_bits(jj,:,:)=HBF_WLLS(K,M,iters,pow_mat,H_noisy,H_nonoise,bit);
    sum_rate_iter_WLLS_2_bits=squeeze(sum_rate_WLLS_2_bits(jj,:,:))
    
    
    
    % -- our method -- %
    
    % noise_cond
    % false: Noiseless channels
    % true: Noisy channels 
    noise_cond=true; 
    
    for noise_idx=1:2
        if noise_idx==1
            noise_cond=false; 
        else
            noise_cond=true; 
        end
            
    % Switch_mat : Optimal switch configuration 
    [sum_rate_our(jj,noise_idx,:), Switch_mat]=Binary_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise,mu,r,divide_rate,noise_cond);
    
    end
    sum_rate_iter_our=squeeze(sum_rate_our(jj,:,:))
end


