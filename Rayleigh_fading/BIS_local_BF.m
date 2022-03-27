function [SINR_max,B_opt]= BIS_local_BF (M,N,L,ini_switch,loc_matrix,All_values,H,P_r,N_P, SINR_max)

for i1=1:size(loc_matrix,1)
    subset=loc_matrix(i1,:);

    for j1=0:N^L-1
        baseStr = dec2base(j1,N,L);
        G1 = baseStr - '0'+1;

        G1=All_values(G1);


        B=ini_switch;
        B(subset)=G1;

        if sum(B==0)~=M
            sol=H*transpose(B);  
            desired_pow=P_r(1)*abs(sol(1).^2);
            Intrf_pow=sum(P_r(2:end)'.*(abs(sol(2:end)).^2));
            n_p=(M-sum(B==0))*N_P;
            SINR=desired_pow/(n_p+Intrf_pow);
        else
            SINR=0;
        end

        if SINR>=SINR_max(i1)
            SINR_max(i1)=SINR;
            B_opt(i1,:)=B;
        end
    end
end