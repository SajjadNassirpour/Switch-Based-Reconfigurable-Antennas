% -- Simplified exhaustive search (SES) with {0,1} -- %
% K :                  Number of users
% M :                  Number of elements in array antenna
% iters :              Number of trials
% pow_mat :            SNR(dB) 
% H_nonoise :          Noiseless channels
% H_noisy :            Noisy channels
% bunch_size:          size of the BF search 
% noisy_cond :         Indicate using of noisy or noiseless channels
% sum_rate_BF:         Sum-rate
% sim_BF_complexity:   Complexity

function [sum_rate_BF, sim_BF_complexity]=Simplified_BF_approach(K,M,bunch_size,iters,pow_mat,H_noisy,H_nonoise,noise_cond)
    
  
    % Average sum-rate
    sum_rate_BF=zeros(1,length(pow_mat));  
        
    H_correct_channel=H_nonoise;

    if noise_cond()
        H_noisy_channel=H_noisy;
    else
        H_noisy_channel=H_nonoise;    
    end


    for pow_idx=1:length(pow_mat)

        P_sym_dB=pow_mat(pow_idx);
        P_sym=10.^(P_sym_dB/10);
        P_r=P_sym*ones(1,K);
        
        % Noise power
        N_P=1; 

        rate_BF=zeros(1,iters);

        for iter=1:iters
            complexity=0;


            correct_channel=squeeze(H_correct_channel(iter,:,:));
            noisy_channel=squeeze(H_noisy_channel(iter,:,:));

            H_BTS=transpose(noisy_channel(1,:));

            H_i=noisy_channel;
            H_i(1,:)=[];
            H_i=transpose(H_i);
            
            % Noisy channel
            H=transpose([H_BTS  H_i]);  

            
            % ------ Brute-Force method ------- %
            
            All_values = [0 1];
            N=length(All_values);
            
            bunch_mat=nchoosek(1:M,bunch_size);
            
            SINR_max=0;
            if M<=50
            r_rand=randi([1 2^M-1],1,1);
            else
                r_rand=randi([1 2^50-1],1,1);
            end
            baseStr = dec2base(r_rand,2,M);
            B_ini = baseStr - '0';
            
            B_ini=zeros(1,M);
            
            for bunch_idx=1:size(bunch_mat,1)
                
                changed_antenna=bunch_mat(bunch_idx,:);
                
                for i1=0:N^bunch_size-1
                    baseStr = dec2base(i1,N,bunch_size);
                    G1 = baseStr - '0'+1;
                    B=B_ini;
                    B(changed_antenna)=All_values(G1);


                    if sum(B==0)~=M
                        complexity=complexity+1;
                        sol=H*transpose(B);  
                        desired_pow=P_r(1)*abs(sol(1).^2);
                        Intrf_pow=sum(P_r(2:end)'.*(abs(sol(2:end)).^2));
                        n_p=(M-sum(B==0))*N_P;
                        SINR=desired_pow/(n_p+Intrf_pow);
                    else
                        SINR=0;
                    end

                    if SINR>SINR_max
                        SINR_max=SINR;
                        B_opt=B;
                    end
                end
            end

            % ---------------------------------------------
            
            
            % Optimal switch configuration

            H_BTS_correct=transpose(correct_channel(1,:));
            H_i_correct=correct_channel;
            H_i_correct(1,:)=[];
            H_i_correct=transpose(H_i_correct);


            % Noiseless channel
            H_correct=transpose([H_BTS_correct  H_i_correct]);

            % SINR
            sol=H_correct*transpose(B_opt);  
            desired_pow=P_r(1)*abs(sol(1).^2);
            Intrf_pow=sum(P_r(2:end)'.*(abs(sol(2:end)).^2));
            n_p=(M-sum(B_opt==0))*N_P;
            SINR_BF=desired_pow/(n_p+Intrf_pow);


            % Rate at each user
            rate_BF(1,iter)=0.5*log2(1+SINR_BF);
            complexity_mat(1,iter)=complexity;
        end
        sim_BF_complexity(1,pow_idx)=mean(complexity_mat);
        sum_rate_BF(1,pow_idx)=K*sum(mean(rate_BF,2));
    end
    

end




      
