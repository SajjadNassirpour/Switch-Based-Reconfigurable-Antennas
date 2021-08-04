% -- Beamnulling + Analog beamforming appraoch (real-world complex channels) -- %
% K :            Number of users
% M :            Number of elements in array antenna
% iters :        Number of trials
% pow_mat :      SNR(dB) 
% H_nonoise :    Noiseless channels
% null_vector :  Beamnulling + Analog beamforming weighting vector 

function [sum_rate]=Beam_null_app_complex(K,M,iters,pow_mat,H_nonoise,null_vector)
    % average sum_rate
    sum_rate=zeros(1,length(pow_mat));  
    % noiseless channel
    H_correct_channel=H_nonoise;  

    for pow_idx=1:length(pow_mat)
        
        SNR=10^(pow_mat(pow_idx)/10); 
        % Rate at each user
        rate=zeros(K,iters);  

        for iter=1:iters
            for BTS=1:K
            correct_channel=squeeze(H_correct_channel(iter,BTS,:,:));
            H_hat_correct=transpose(correct_channel(BTS,:));  % Desired channel
            H_i_correct=correct_channel;
            H_i_correct(BTS,:)=[];
            H_i_correct=transpose(H_i_correct);  % Interference channel

            % weitghing vector
            W_1=transpose(squeeze(null_vector(iter,BTS,:)));  
            
            % Desired signal
            Desired_sig=W_1*H_hat_correct;
            % Interference signal
            intrf_sig=W_1*H_i_correct;

            % Generate noise at receiver
            a = 1/sqrt(SNR); 
            n_R=(a.*randn(1,M))./sqrt(2);
            n_I=(a.*randn(1,M))./sqrt(2);
            n=n_R+1i*n_I;
            
            % Desired power
            Desired_pow=abs(Desired_sig)^2;
            % Interference power
            Intrf_power=sum(abs(intrf_sig).^2); 
            % Noise power
            noise_pow=sum(abs(W_1.*n).^2);  
            
             % SINR
            SINR=Desired_pow/(noise_pow+Intrf_power); 
            
            % Rate at each user
            rate(BTS,iter)=log2(1+SINR);   
            end
        end
        % Sum-rate
        sum_rate(1,pow_idx)=sum(mean(rate,2));
    end
    
end









