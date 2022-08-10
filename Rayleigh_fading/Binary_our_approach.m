% -- Our optimization appraoch based on sigmoid filled function -- %

% K :              Number of users
% M :              Number of elements in array antenna
% iters :          Number of trials
% pow_mat :        SNR(dB) 
% H_nonoise :      Noiseless channels
% H_noisy :        Noisy channels
% r,divide_rate:   Optimization parameters
% noisy_cond :     Indicate using of noisy or noiseless channels
% our_final_rate:  Sum-rate
% BF_our_ratio :   The Complexity (# of evaluations) of our approach
% TDMA_rate :      Sum-rate with TDMA approach

function [our_final_rate, Complexity_ratio,TDMA_rate]=Binary_our_approach(K,M,iters,H_noisy,H_nonoise,pow_mat,mu,r,divide_rate,noisy_cond,initial_rand_num)

    pow_mat=10.^(pow_mat/10);
    P_r=pow_mat*ones(1,K); 
    
    % Noise power
    N_P=1; 

    % Optimization Parameters
    input_mu=mu;
    input_r=r;

    % SINR and rate vectors
    SINR_GS=zeros(1,iters);
    indvidual_rate=zeros(1,iters);
    
    % Number of evaluations
    Eval=zeros(1,iters);
    
    for iter=1:iters
        
        eval1=1;
        
        % The channels
        if (noisy_cond)
           H=squeeze(H_noisy(iter,:,:));
        else
           H=squeeze(H_nonoise(iter,:,:)); 
        end
        min_f_total=0;
        for internal_iter=1:initial_rand_num
            % Generate random initial switch configuration
            All_values=[0 1];
            N=length(All_values);

            ar = randi([1,N],1,M);
            G1_ini = All_values(ar);

            %===== Filled-function approach =====%
            mu=input_mu;
            r=input_r;

            r0=r;
            resume=1;

            step=0;
            ii=0;
            kk=0;
            S_star_mat=[];
            S_final=[];
            min_f=0;
            while resume==1
                step=step+1;
                jump=0;   
                [S_star,eval]=LS_algo_diff_pow(G1_ini,H,P_r,N_P,1,r);
                eval1=eval1+eval;
                S_star_mat(kk+1,:)=S_star;

                if kk==0
                    if obj_func_diff_pow(S_star_mat(kk+1,:),H,P_r,N_P)>= min_f
                    min_f=obj_func_diff_pow(S_star_mat(kk+1,:),H,P_r,N_P);
                    S_final=S_star_mat(kk+1,:);
                    end
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
            if min_f<min_f_total
                min_f_total=min_f;
                S_final_total=S_final;
            end
            
        end

        % Final evaluation based on optimal switch configuration
        B_opt_GS=S_final_total;
        H_correct=squeeze(H_nonoise(iter,:,:));
        SINR_GS(iter)=-obj_func_diff_pow2(B_opt_GS,H_correct,P_r,N_P);
        indvidual_rate(iter)=0.5*log2(1+SINR_GS(iter));
        
        % Number of evaluations
        Eval(iter)=eval1;
    end
    
    % Average results
    Complexity_ratio=mean(Eval);
    TDMA_rate=(1/2)*log2(1+P_r(1));
    our_final_rate=K*mean(indvidual_rate);
end


      
