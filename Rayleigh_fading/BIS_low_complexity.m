% -- BIS genetic algorithm when s_m \in {0,1} -- %
% K :              Number of users
% M :              Number of elements in array antenna
% iters :          Number of trials
% pow_mat :        SNR(dB) 
% H_nonoise :      Noiseless channels
% H_noisy :        Noisy channels
% Ns:              Number of survived parents
% L:               Number of elements for BF in subsets
% noisy_cond :     Indicate using of noisy or noiseless channels

% final_sum_rate: Sum-rate
% Switch_mat :     The matrix of the optimal switch confgurations

function [final_sum_rate, Switch_mat]=BIS_low_complexity(K,M,iters,pow_mat,H_noisy,H_nonoise,Ns,L,noise_cond)
    
    Switch_mat=zeros(iters,length(pow_mat),1,M);
    % Average sum-rate
    final_sum_rate=zeros(1,length(pow_mat));  
        

    

    for pow_idx=1:length(pow_mat)

        P_sym_dB=pow_mat(pow_idx);
        P_sym=10.^(P_sym_dB/10);
        P_r=P_sym*ones(1,K);
        
        % Noise power
        N_P=1; 

        rate_BIS=zeros(K,iters);

        for iter=1:iters
            
            % The channels
            if (noise_cond)
               H=squeeze(H_noisy(iter,:,:));
            else
               H=squeeze(H_nonoise(iter,:,:)); 
            end

            
            % ------ BIS low-complexity method ------- %
            
            All_values = [0 1];
            N=length(All_values);
            loc_matrix=nchoosek(1:M,L);
            ini_switch=zeros(1,M);
            SINR_max=zeros(size(loc_matrix,1),1);
            step1=1;
            
            [SINR_max,B_opt]= BIS_local_BF (M,N,L,ini_switch,loc_matrix,All_values,H,P_r,N_P, SINR_max);
            
            fix_loc_first=[];    
            [sort_SINR_max, ia, ic] = unique(SINR_max);
           

            B_0_new=B_opt(ia(end-Ns+1:end),:);
            SINR_max_new=SINR_max(ia(end-Ns+1:end));

            fix_loc_first= loc_matrix(ia(end-Ns+1:end),:);
                       
            while L*step1<M
                step1=step1+1;
                B_opt_total=[];
                SINR_max_total=[];
                loc_matrix_total=[];
                for i2=1:Ns
                    
                    
                    M_remained=1:M;
                    
                    M_remained(fix_loc_first(i2,:))=[];
                    
                    loc_matrix=nchoosek(M_remained,L);
                    ini_switch=B_0_new(i2,:);
                    SINR_max=SINR_max_new(i2)*ones(size(loc_matrix,1),1);
                    
                    [SINR_max,B_opt]= BIS_local_BF (M,N,L,ini_switch,loc_matrix,All_values,H,P_r,N_P, SINR_max);

                    B_opt_total=[B_opt_total; B_opt];
                    SINR_max_total=[SINR_max_total; SINR_max];
                    loc_matrix_total=[loc_matrix_total; [repmat(fix_loc_first(i2,:),size(loc_matrix,1),1) loc_matrix]];
                    
                end
                
                [sort_SINR_max, ia, ic] = unique(SINR_max_total);
                step1;
                
                if length(ia)==1
                    
                    ia=repmat(ia,Ns,1);
                    
                end
                B_0_new=B_opt_total(ia(end-Ns+1:end),:);
                SINR_max_new=SINR_max_total(ia(end-Ns+1:end));

                fix_loc_first= loc_matrix_total(ia(end-Ns+1:end),:);
            end
                
            
            [sort_SINR_max, ia, ic] = unique(SINR_max_total);
            B_opt=B_opt_total(ia(end),:);


            % Optimal switch configuration
            B_opt_BIS=B_opt;  

            Switch_mat(iter,pow_idx,1,:)=B_opt_BIS;
            H_correct=squeeze(H_nonoise(iter,:,:));
            % SINR
            sol=H_correct*transpose(B_opt_BIS);  
            desired_pow=P_r(1)*abs(sol(1).^2);
            Intrf_pow=sum(P_r(2:end)'.*(abs(sol(2:end)).^2));
            n_p=(M-sum(B_opt_BIS==0))*N_P;
            SINR_BIS=desired_pow/(n_p+Intrf_pow);


            % Rate at each user
            rate_BIS(1,iter)=0.5*log2(1+SINR_BIS);

        end
    % Average sum-rate
    final_sum_rate(1,pow_idx)=K*sum(mean(rate_BIS,2));

    end
    
 
end




      
