% lattice-based method

function FRF= Alt_SIC(Fopt,cons,NRF)
    Fopt1=(Fopt-mean(Fopt))/std(Fopt);
    for rf=1:NRF
        for i=1:length(Fopt)
            if cons==0
                FRF(i)=ClimbOne(Fopt(i),cons);
            else
                FRF(i)=ClimbOne(Fopt1(i),cons);
            end
        end
    end
    ss=1;
end


function s_out=ClimbOne(s_in,cons)

    if cons==0
        s_out=s_in;
    else
        len=size(cons,2);

        s_com2=ones(1,len)*s_in;

        s_com=abs(cons-s_com2);

        ind=find(s_com==min(s_com));

        if ~isempty(ind)
            s_out=cons(ind(1));
        else
            s_out=1;
        end
    end
end


