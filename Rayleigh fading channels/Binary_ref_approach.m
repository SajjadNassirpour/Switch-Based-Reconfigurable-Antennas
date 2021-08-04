% -- Optimization appraoch based on the following filled function -- %

% Y. H. Gu and Z. Y. Wu, “A new filled function method for nonlinear 
% integer programming problem,”Applied mathematics and computation,
% vol. 173, no. 2, pp. 938–950, 2006.

% K :              Number of users
% M :              Number of elements in array antenna
% iters :          Number of trials
% pow_mat :        SNR(dB) 
% H_nonoise :      Noiseless channels
% r,divide_rate:   Optimization parameters

% our_final_rate:  Sum-rate
% BF_ref_ratio :   The tatio of (# of evaluations of Brute Force/ the ref approach)

function [ref_final_rate, BF_ref_ratio]=Binary_ref_approach(K,M,iters,H_nonoise,pow_mat,mu,r,divide_rate)

    pow_mat=10.^(pow_mat/10);
    P_r=pow_mat*ones(1,K);

    % Noise power
    N_P=1; 

    % Optimization Parameters
    input_mu=mu;
    input_r=r;

    % SINR matrix 
    SINR_LN_rand=zeros(1,iters);
    
    % Number of evaluations
    Eval=zeros(1,iters);

    for iter=1:iters
        
        iter
        eval1=1;
        
        % The channels
        H=squeeze(H_nonoise(iter,:,:));
        
        % Generate random initial switch configuration        
        if M<=50
        r_rand=randi([1 2^M-1],1,1);
        else
            r_rand=randi([1 2^50-1],1,1);
        end
        baseStr = dec2base(r_rand,2,M);
        G1_ini = baseStr - '0';
        
        %===== Filled-function approach =====%
        mu=input_mu;
        r=input_r;
        
        r0=r;
        resume=1;
        
        step=0;
        ii=0;
        kk=0;
        S_star_mat=[];
        while resume==1
            step=step+1;
            jump=0;   
            [S_star,eval]=LS_algo_ref_diff_pow(G1_ini,H,P_r,N_P,1,r);
            eval1=eval1+eval;
            S_star_mat(kk+1,:)=S_star;
            
            if kk==0
                min_f=obj_func_diff_pow(S_star_mat(kk+1,:),H,P_r,N_P);
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
                [S_bar,eval]=LS_algo_ref_diff_pow(S_star,H,P_r,N_P,2,r);
                eval1=eval1+eval;
                G1_ini = S_bar;
                kk=kk+1;
            end

            if ii>=M+1 
                if r>=mu
                    r=r/divide_rate;
                    ii=1;
                    S_star=S_star_mat(kk+1,:);
                    [S_bar,eval]=LS_algo_ref_diff_pow(S_star,H,P_r,N_P,2,r);
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
                [S_bar,eval]=LS_algo_ref_diff_pow(S_star_j0,H,P_r,N_P,2,r);
                eval1=eval1+eval;
                G1_ini = S_bar;
                kk=kk+1;
            end
        end

        B_opt_LN_rand=S_final;
        score_LN_rand=-obj_func_diff_pow(B_opt_LN_rand,H,P_r,N_P);

        SINR_LN_rand(iter)=score_LN_rand;
        Eval(iter)=eval1;
    end
    
    % Average results
    LNN_rand=mean(SINR_LN_rand);
    Avg_Eval=mean(Eval);
    indvidual_rate=(1/2)*log2(1+LNN_rand);

    BF_ref_ratio=(2^M)/mean(Avg_Eval);
    ref_final_rate=K*indvidual_rate;
end


      
