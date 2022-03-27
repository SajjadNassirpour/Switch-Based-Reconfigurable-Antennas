% -- Lattice-base beamforming (real-world complex channels) -- %


function [sum_rate]=HBF_Lattice(K,M,iters,pow_mat,H_noisy,H_nonoise,bit,ee)

N_RF=1;

n0=bit; n=2^n0;% n0 is the number of bits
cons=exp(1i*[0:2*pi/n:2*pi-2*pi/n]);% constellation without '0'
cons=[cons,0];%constellation with '0'

if bit==0
    cons=0;
end
if bit==1
    cons=[0 1];
end

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
            rate=zeros(K,iters);
            for iter=1:iters
                for BTS=1:K
                correct_channel=squeeze(H_correct_channel(iter,BTS,:,:));
                noisy_channel=squeeze(H_noisy_channel(iter,BTS,:,:));

                % Correct channels
                H_hat_correct=transpose(correct_channel(BTS,:));
                H_i_correct=correct_channel;
                H_i_correct(BTS,:)=[];
                H_i_correct=transpose(H_i_correct);

                % Noisy channels
                H_hat=transpose(noisy_channel(BTS,:));
                H_i=noisy_channel;
                H_i(BTS,:)=[];
                H_i=transpose(H_i);
                R=zeros(M);
                for intrf_idx=1:K-1
                    H_i_1=H_i(:,intrf_idx);
                    R=R+H_i_1*(H_i_1');
                end
                
                R=R+ee*eye(M)/SNR;
                
                if ee==0
                    inv_mat=pinv(H_hat*(H_hat')+R);
                else
                    inv_mat=inv(H_hat*(H_hat')+R);
                end
                W_1=(H_hat')*inv_mat;
                
                
                W_1=W_1/norm(W_1,'fro');
                
                Fopt=transpose(W_1);
               
                W_1 =Alt_SIC(Fopt,cons,N_RF);

                

                

                W_1=W_1/norm(W_1,'fro');

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
                Desired_pow=sum(abs(Desired_sig).^2);
                % Interference power
                Intrf_pow=sum(abs(intrf_sig).^2);
                % Noise power
                noise_pow=sum(abs(W_1.*n).^2);
                
                % SINR
                SINR=Desired_pow/(noise_pow+Intrf_pow);
                
                % Rate at each user
                rate(BTS,iter)=log2(1+SINR);
                end
            end
            % Sum-rate
            sum_rate(noise_cond+1,pow_idx)=sum(mean(rate,2));
        end
    end
end









