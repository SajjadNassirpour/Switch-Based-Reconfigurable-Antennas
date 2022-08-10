% -- Successive refinemnet(SR) method with {0,1} -- %
% K :              Number of users
% M :              Number of elements in array antenna
% iters :          Number of trials
% pow_mat :        SNR(dB) 
% H_nonoise :      Noiseless channels
% H_noisy :        Noisy channels
% threshold:       Stopping condition
% noisy_cond :     Indicate using of noisy or noiseless channels
% sum_rate_SR:     Sum-rate
% complexity_SR:   Complexity

function [sum_rate_SR, complexity_SR]=SR_approach(K,M,iters,pow_mat,H_noisy,H_nonoise,noise_cond,threshold,rand_ini_num)
    
  
    % Average sum-rate
    sum_rate_SR=zeros(1,length(pow_mat));  
    complexity_SR=zeros(1,length(pow_mat)); 
        
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

        rate_SR=zeros(1,iters);

        for iter=1:iters

            correct_channel=squeeze(H_correct_channel(iter,:,:));
            noisy_channel=squeeze(H_noisy_channel(iter,:,:));

            H_BTS=transpose(noisy_channel(1,:));

            H_i=noisy_channel;
            H_i(1,:)=[];
            H_i=transpose(H_i);
            
            % Noisy channel
            H=transpose([H_BTS  H_i]);  

            
            % ------ SR method ------- %
            
            All_values = [0 1];
            N=length(All_values);
            
            eval_mat=0;
            SINR_max=0;
            
            for rand_idx=1:rand_ini_num

                
                ar = randi([1,N],1,M);
                B_ini = All_values(ar);
                 
            
                SINR_ini=0;
                SINR_thre=inf;
                
                while sum(abs(SINR_ini-SINR_thre))>=threshold

                    SINR_thre=SINR_ini;
                    for i1=1:M
                        for item=1:N
                            if B_ini(1,i1)~=All_values(item) 
                                eval_mat=eval_mat+1;
                                B=B_ini;
                                B(1,i1)=All_values(item);
                                if sum(B==0)~=M
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
                                    SINR_opt=SINR;
                                    B_opt=B;
                                end
                            end
                        end
                    end
                    if SINR_opt==SINR_ini
                    else
                        SINR_ini=SINR_opt;
                        B_ini=B_opt;
                    end
                end

            end
            
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
            SINR_SR=desired_pow/(n_p+Intrf_pow);


            % Rate at each user
            rate_SR(1,iter)=0.5*log2(1+SINR_SR);
        end
        complexity_SR(1,pow_idx)=mean(eval_mat); 
        sum_rate_SR(1,pow_idx)=K*sum(mean(rate_SR,2));
         
    end

end