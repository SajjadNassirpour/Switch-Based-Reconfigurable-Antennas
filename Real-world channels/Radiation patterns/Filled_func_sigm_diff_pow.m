% -- The sigmoid filled function -- %
% S:              Switch configuration
% S_star:         Optimal switch configuration
% H:              Channel coefficients
% t_p:            Transmitted power
% n_p:            Noise power
% t:              f(S)-f(S_star)
% r:              Optimization parameter
% F_r:            Filled finction value

function F_r=Filled_func_sigm_diff_pow(S_star,S,r,H,t_p,n_p)

    f_S=obj_func_diff_pow(S,H,t_p,n_p);
    f_S_star=obj_func_diff_pow(S_star,H,t_p,n_p);
    t=f_S-f_S_star;
    
    %----g_r(t)-----%
    if t>=0
        g=1;
    elseif t<=-r
        g=t+r;
    else
        g=1/(1+exp((-6/r)*(t+r/2)));
    end
    
    %----f_r(t)-----%
    if t>=0
        f=1;
    elseif t<=-r
        f=t+r;
    else
        f=1/(1+exp((-6/r)*(t+r/2)));
    end
    
    %----Filled function-----%
    
    if t<=-r
        cond=0;
    else
        cond=1;
    end
    
    F_r=(1/(cond*(norm(S-S_star)^2)+1))*g+f;
end