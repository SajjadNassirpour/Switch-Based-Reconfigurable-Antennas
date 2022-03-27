% -- Objective function of the optimzaion problem -- %
% S:              Switch configuration
% H:              Channel coefficients
% t_p:            Transmitted power
% n_p:            Noise power
% snir_eval:      SINR value

function snir_eval=obj_func_diff_pow(S,H,t_p,n_p)
    M=length(S);
    if sum(S==0)~=M
        S1=S/norm(S);
        B=S1';
        sol=H*B;  
        desired_pow=t_p(1)*abs(sol(1).^2);
        Intrf_pow=sum(t_p(2:end)'.*(abs(sol(2:end)).^2));
        n_p=sum(S==1)*n_p/(norm(S)^2);
        snir_eval=-(desired_pow/(n_p+Intrf_pow));
    else
        snir_eval=0;
    end
end