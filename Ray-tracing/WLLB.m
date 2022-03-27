% Calculate WLLS weighting vectors

function F =WLLB(Frf_temp,Wrf_temp,H,Nt,Nr,Nrf,b,iter)

if nargin <=7
    iter=10;
end

if b==0
    F=transpose(H);
else
    
    omiga=exp(1j*2*pi/(2^b));
    f_ind=[0:1:2^b-1]';
    F_set=omiga.^(f_ind); 
    
    if b==1
        F_set=[0; 1];
    else  
        F_set=[0; F_set];
    end

    Q=H;
    [U S V]=svd(H);
    U=U(:,1:Nrf);
    V=V(:,1:Nrf);
    S=S(1:Nrf,1:Nrf);

    V=transpose(H);

    ff=zeros(Nt,Nrf);
    ww=zeros(Nr,Nrf);

    Frf=zeros(Nt,Nrf);
    Wrf=zeros(Nr,Nrf);

    Mat_f=zeros(Nt,size(F_set,1));
    Mat_w=zeros(Nr,size(F_set,1));


    sigma=1;
    ncount=0;
    while sigma >0.01 && ncount<=iter
        ncount=ncount+1;
        for k=1:Nrf
            ff(:,k)=[];
            ww(:,k)=[];
            Q=U*inv(0.01*eye(Nrf)+V'*ff*ww'*U)*V'*S;      

            for i=1:Nt
                for ii=1:size(F_set,1)
                    Mat_f(:,ii) = Frf_temp(:,k);
                    Mat_f(i,ii)=F_set(ii);
                end


                Product = Wrf_temp(:,k)'*Q*Mat_f;

                [~, position] = max(Product);
                Frf_temp(:,k) = Mat_f(:,position);
            end
            ff=Frf_temp;
            ww=Wrf_temp;
        end
        sigma1=norm(Frf-Frf_temp,'fro');
        sigma2=norm(Wrf-Wrf_temp,'fro');
        sigma=max(sigma1,sigma2);
        Frf=Frf_temp;
        Wrf=Wrf_temp;
    end

    F=Frf;
end






