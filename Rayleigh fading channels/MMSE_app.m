% -- MMSE appraoch -- %
% K :         Number of users
% M :         Number of elements in array antenna
% iters :     Number of trials
% pow_mat :   SNR(dB) 
% H_nonoise : Noiseless channels
% H_noisy :   Noisy channels

function [sum_rate]=MMSE_app(K,M,iters,pow_mat,H_noisy,H_nonoise)

    sum_rate=zeros(2,length(pow_mat));
    for noise_cond=0:1
        % noise_cond = 1 means noisy channels are used
        if noise_cond==0
            H_correct_channel=H_nonoise;
            H_noisy_channel=H_nonoise;
        else
            H_correct_channel=H_nonoise;
            H_noisy_channel=H_noisy;
        end

        for pow_idx=1:length(pow_mat)
            SNR=10^(pow_mat(pow_idx)/10);
            rate=zeros(1,iters);
            for iter=1:iters

                correct_channel=squeeze(H_correct_channel(iter,:,:));
                noisy_channel=squeeze(H_noisy_channel(iter,:,:));
                
                % Correct channels
                H_hat_correct=transpose(correct_channel(1,:));
                H_i_correct=transpose(correct_channel(2:end,:));
                
                % Noisy channels
                H_hat=transpose(noisy_channel(1,:));
                H_i=transpose(noisy_channel(2:end,:)); % Interference channels

                R=H_i*(H_i')+eye(M)/SNR;
                inv_mat=inv(H_hat*(H_hat')+R);
                W_1=(H_hat')*inv_mat;
                % Desired signal
                Desired_sig=W_1*H_hat_correct;
                
                % Interference signal
                intrf_sig=W_1*H_i_correct;
                
                % Generate noise at receiver
                a = 1/sqrt(SNR); 
                n=a.*randn(1,M);
                
                % Desired power
                Desired_pow=abs(Desired_sig)^2;
                % Interference power
                Intrf_pow=sum(abs(intrf_sig).^2);
                % Noise power
                noise_pow=sum(abs(W_1.*n).^2);
                
                % SINR
                SINR=Desired_pow/(noise_pow+Intrf_pow);
                rate(iter)=1/2*log2(1+SINR);
            end
            % Sum-rate
            sum_rate(noise_cond+1,pow_idx)=K*mean(rate);
        end
    end
end









