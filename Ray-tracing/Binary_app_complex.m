% -- Our optimization approach based on filled function (real-world complex channel) -- %
% K :              Number of users
% M :              Number of elements in array antenna
% iters :          Number of trials
% pow_mat :        SNR(dB) 
% H_nonoise :      Noiseless channels
% H_noisy :        Noisy channels
% r,divide_rate:   Optimization parameters
% noisy_cond :     Indicate using of noisy or noiseless channels

% final_sum_rate:  Sum-rate
% Switch_mat :     The matrix of the optimal switch confgurations

function [final_sum_rate, Switch_mat]=Binary_app_complex(K,M,iters,pow_mat,H_noisy,H_nonoise,mu,r,divide_rate,noise_cond)
    
    Switch_mat=zeros(iters,length(pow_mat),K,M);
    
    initial_r=r;
    % Average sum-rate
    final_sum_rate=zeros(1,length(pow_mat));  
        
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

        rate_GS=zeros(K,iters);

        for iter=1:iters

            for BTS=1:K

            eval1=1;
            correct_channel=squeeze(H_correct_channel(iter,BTS,:,:));
            noisy_channel=squeeze(H_noisy_channel(iter,BTS,:,:));

            H_BTS=transpose(noisy_channel(BTS,:));

            H_i=noisy_channel;
            H_i(BTS,:)=[];
            H_i=transpose(H_i);
            
            % Noisy channel
            H=transpose([H_BTS  H_i]);  

            % Generate random initial switch configuration
            if M<=50
            r_rand=randi([1 2^M-1],1,1);
            else
                r_rand=randi([1 2^50-1],1,1);
            end
            baseStr = dec2base(r_rand,2,M);
            G1_ini = baseStr - '0';
            
            %===== Filled-function approach =====%
            r=initial_r;
            r0=r;
            resume=1;
            
            step=0;
            ii=0;
            kk=0;
            S_star_mat=[];
            while resume==1
                step=step+1;
                jump=0;   
                % LS_algo_diff_pow ---> local and global search
                [S_star,eval]=LS_algo_diff_pow(G1_ini,H,P_r,N_P,1,r);
                eval1=eval1+eval;
                S_star_mat(kk+1,:)=S_star;
                
                % obj_func_diff_pow ---> SINR as objective function
                if kk==0
                    min_f=obj_func_diff_pow(S_star_mat(kk+1,:),H,P_r,N_P);
                    % S_final --> optimal switch configuration
                    S_final=S_star_mat(kk+1,:);
                end
                
                if kk>0
                   if obj_func_diff_pow(S_star_mat(kk+1,:),H,P_r,N_P)>= min_f
                       jump=1;
                   else
                       min_f=obj_func_diff_pow(S_star_mat(kk+1,:),H,P_r,N_P);
                       S_final=S_star_mat(kk+1,:);
                       r=r0;
                   end
                end

                if jump==0
                    ii=1;
                    S_star=S_star_mat(kk+1,:);
                    [S_bar,eval]=LS_algo_diff_pow(S_star,H,P_r,N_P,2,r);
                    eval1=eval1+eval;
                    G1_ini = S_bar;
                    kk=kk+1;
                end

                if ii>=M+1 
                    if r>=mu
                        r=r/divide_rate;
                        ii=1;
                        kk=kk-1;
                        S_star=S_star_mat(kk+1,:);
                        [S_bar,eval]=LS_algo_diff_pow(S_star,H,P_r,N_P,2,r);
                        eval1=eval1+eval;
                        G1_ini = S_bar;
                        kk=kk+1;
                        jump=0;
                    else
                        resume=0;
                    end
                end    

                if ii<M+1 && jump==1
                    ii=ii+1;
                    kk=kk-1;
                    S_star=S_star_mat(kk+1,:);
                    G1_x=zeros(1,M);
                    G1_x(1,ii-1)=1;
                    S_star_j0=double(xor(S_star,G1_x));
                    [S_bar,eval]=LS_algo_diff_pow(S_star_j0,H,P_r,N_P,2,r);
                    eval1=eval1+eval;
                    G1_ini = S_bar;
                    kk=kk+1;
                end

            end
            % Optimal switch configuration
            B_opt_GS=S_final;  

            Switch_mat(iter,pow_idx,BTS,:)=B_opt_GS;

            H_BTS_correct=transpose(correct_channel(BTS,:));
            H_i_correct=correct_channel;
            H_i_correct(BTS,:)=[];
            H_i_correct=transpose(H_i_correct);


            % Noiseless channel
            H_correct=transpose([H_BTS_correct  H_i_correct]);

            % SINR
            SINR_GS=-obj_func_diff_pow2(B_opt_GS,H_correct,P_r,N_P);

            % Rate at each user
            rate_GS(BTS,iter)=log2(1+SINR_GS);
            end
        end
    % Average sum-rate
    final_sum_rate(1,pow_idx)=sum(mean(rate_GS,2));
    end 
end




      
