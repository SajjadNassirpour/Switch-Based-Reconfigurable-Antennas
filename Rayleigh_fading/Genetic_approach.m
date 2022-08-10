% -- BF method when with {0,1} -- %
% K :              Number of users
% M :              Number of elements in array antenna
% iters :          Number of trials
% pow_mat :        SNR(dB) 
% H_nonoise :      Noiseless channels
% H_noisy :        Noisy channels
% noisy_cond :     Indicate using of noisy or noiseless channels

% sum_rate_gen:    Sum-rate
% GA_complexity :  Complexity

function [sum_rate_gen,GA_complexity]=Genetic_approach(K,M,pop_size,num_point_cross,mut_prob,iters,pow_mat,H_noisy,H_nonoise,noise_cond)
    
    % Number of iterations of the GA method until convergence
    GA_iters=100;
  
    % Average sum-rate
    sum_rate_gen=zeros(1,length(pow_mat));  
        
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
            
            SINR_max=0;
            
            for pop_idx=1:pop_size
                G1_rand=randi(1:N,1,M);
                B=All_values(G1_rand);
                B_mat(pop_idx,:)=B;
            end
            for pop_idx=1:pop_size
                B=B_mat(pop_idx,:);

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
                
                SINR_mat(pop_idx)=SINR;
            end
            if max(SINR_mat)>SINR_max
                [SINR_max,SINR_max_loc]=max(SINR_mat);
                B_opt=B_mat(SINR_max_loc,:);
            end
            stop=0;
            for itt=1:GA_iters
            
                Actual_count_mat_new=zeros(1,pop_size);
                for pop_idx=1:pop_size
                    avg_SINR=mean(SINR_mat);
                    Prob_mat(pop_idx)=SINR_mat(pop_idx)/sum(SINR_mat);
                    Acum_prob_mat(pop_idx)=sum(Prob_mat(1:pop_idx));

                    Expected_count_mat(pop_idx)=SINR_mat(pop_idx)/avg_SINR;
                    Actual_count_mat(pop_idx)=round(Expected_count_mat(pop_idx));

                    rand_num=rand;
                    diff0=Acum_prob_mat-rand_num;
                    diff1=100*ones(1,pop_size);
                    diff1(diff0>=0)=diff0(diff0>=0);
                    [min_val,min_loc]=min(diff1);
                    Actual_count_mat_new(min_loc)=Actual_count_mat_new(min_loc)+1;

                    B_mat_ini_cross(pop_idx,:)=B_mat(min_loc,:);
                end
                
                

                for cross_idx=1:pop_size/2
                    rand_int_num=sort(randi([1 M-1],1,num_point_cross));
                    rand_int_num=[rand_int_num M];
                    
                    B_mat_cross(2*(cross_idx-1)+1:2*cross_idx,:)=B_mat_ini_cross(2*(cross_idx-1)+1:2*cross_idx,:);
                
                    for num_p_cross_idx=1:num_point_cross
                        if mod(num_p_cross_idx,2)==0
                            B_mat_cross(2*(cross_idx-1)+1,rand_int_num(num_p_cross_idx)+1:rand_int_num(num_p_cross_idx+1))=B_mat_ini_cross(2*cross_idx,rand_int_num(num_p_cross_idx)+1:rand_int_num(num_p_cross_idx+1));
                        end
                    end
                end
                
                
                rand_mat=rand(pop_size,M);
                B_mat_mut=B_mat_cross;
                for pop_idx=1:pop_size
                    for m_idx=1:M
                        if rand_mat(pop_idx,m_idx)<mut_prob
                           B_mat_mut(pop_idx,m_idx)=1-B_mat_mut(pop_idx,m_idx);
                        end
                    end
                end

                for pop_idx=1:pop_size
                    B=B_mat_mut(pop_idx,:);

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

                    SINR_mat(pop_idx)=SINR;
                end
                
                for pop_idx=1:pop_size
                    B=B_mat_cross(pop_idx,:);

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

                    SINR_mat_cross(pop_idx)=SINR;
                end
                
                
                if max(SINR_mat)>SINR_max
                    [SINR_max,SINR_max_loc]=max(SINR_mat);
                    B_opt=B_mat_mut(SINR_max_loc,:);
                    B_mat=B_mat_mut;
                end
                if max(SINR_mat_cross)>SINR_max
                    [SINR_max,SINR_max_loc]=max(SINR_mat_cross);
                    B_opt=B_mat_cross(SINR_max_loc,:);
                    B_mat=B_mat_cross;
                end
                if max(SINR_mat)<=SINR_max && max(SINR_mat_cross)<=SINR_max 
                    stop=1;
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
        GA_complexity(1,pow_idx)=mean(complexity_mat);
        sum_rate_gen(1,pow_idx)=K*sum(mean(rate_BF,2));
    end
    

end




      
