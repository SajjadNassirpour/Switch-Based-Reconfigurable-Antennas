% -- Local search algorithm based on the following paper-- %

% Y. H. Gu and Z. Y. Wu, “A new filled function method for nonlinear 
% integer programming problem,”Applied mathematics and computation,
% vol. 173, no. 2, pp. 938–950, 2006.

% ini_string:   Initial switch for local search
% H:            Channel coefficients
% t_p:          Transmitted power
% n_p:          Noise power
% LM_string:    Local minimizer switch
% eval :        Number of evaluations

function [LM_string,eval]=LS_algo_ref2_diff_pow(ini_string,H,t_p,n_p,mode,r)
    eval=0;    
    M=length(ini_string);
    S_star=ini_string;
    
    % mode=1: Local search based on SINR
    % mode=2: Local search based on Filled function
    if mode==1
        score_max=obj_func_diff_pow(ini_string,H,t_p,n_p);
    else
        score_max=Filled_func2_diff_pow(S_star,S_star,r,H,t_p,n_p);
    end
    
    % First round of local searching process
    G1_ini = ini_string;
    LM_string=ini_string;
    for i1=1:M
        G1_x=zeros(1,M);
        G1_x(1,i1)=1;
        G1=double(xor(G1_ini,G1_x));
        if sum(G1==0)~=M
            eval=eval+1;
            if mode==1
                score=obj_func_diff_pow(G1,H,t_p,n_p);
            else
                score=Filled_func2_diff_pow(S_star,G1,r,H,t_p,n_p);
            end
            if score<=score_max 
                score_max=score;
                LM_string=G1;
            end
        end
    end
    
    % Run the while loop until G1_ini will be the local optimum
    ok=0;
    while ok==0 
        G1_ini = LM_string;
        for i1=1:M
            G1_x=zeros(1,M);
            G1_x(1,i1)=1;
            G1=double(xor(G1_ini,G1_x));
            if sum(G1==0)~=M
                eval=eval+1;
                if mode==1
                    score=obj_func_diff_pow(G1,H,t_p,n_p);
                else
                    score=Filled_func2_diff_pow(S_star,G1,r,H,t_p,n_p);
                end
                if score<=score_max 
                    score_max=score;
                    LM_string=G1;
                end
            end
        end
        if sum(abs(LM_string-G1_ini)==0)==M 
            ok=1;
        end
    end
end




