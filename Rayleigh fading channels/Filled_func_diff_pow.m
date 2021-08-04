% -- The filled function of the following paper-- %

% Y. H. Gu and Z. Y. Wu, “A new filled function method for nonlinear 
% integer programming problem,”Applied mathematics and computation,
% vol. 173, no. 2, pp. 938–950, 2006.

% S:              Switch configuration
% S_star:         Optimal switch configuration
% H:              Channel coefficients
% t_p:            Transmitted power
% n_p:            Noise power
% t:              f(S)-f(S_star)
% r:              Optimization parameter
% F_r:            Filled finction value

function F_r=Filled_func_diff_pow(S_star,S,r,H,t_p,n_p)

    f_S=obj_func_diff_pow(S,H,t_p,n_p);
    f_S_star=obj_func_diff_pow(S_star,H,t_p,n_p);
    t=f_S-f_S_star;
    
    %----g_r(t)-----%
    if t>0
        g=1;
    elseif t<=-r
        g=0;
    else
        g=((-2)/(r^3))*t^3+((-3)/(r^2))*t^2+1;
    end
    
    %----f_r(t)-----%
    if t>0
        f=1;
    elseif t<=-r
        f=t+r;
    else
        f=((r-2)/(r^3))*t^3+((r-3)/(r^2))*t^2+1;
    end
    
    %----Filled function-----%
    
    F_r=(1/(norm(S-S_star)^2+1))*g+f;
end