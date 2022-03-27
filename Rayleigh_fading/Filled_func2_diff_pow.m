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

function F_r=Filled_func2_diff_pow(S_star,S,r,H,t_p,n_p)

    f_S=obj_func_diff_pow(S,H,t_p,n_p);
    f_S_star=obj_func_diff_pow(S_star,H,t_p,n_p);
    t=f_S-f_S_star;
    
    F_r=-norm(S_star-S)+r*(min(0,t)^3);
    
    
    
end