function [expI, CoM] = expand(I,N) 
%%% Expand the array by even power of 2
    L = length(I);
    expI = I;
    for n=1:N
        if n>1
            L = length(expI);
            I = expI;
        end
        L=L*2;
        expI = zeros(1,L);
        k = 1;
        for j=1:L-1
            if mod(j,2)==1 % odd
                expI(j)=I(k);
                k=k+1;
                if k>L
                    k=k-1;
                end
            else % even
                expI(j)=(I(k-1)+I(k))/2;
            end
        end
    end
    
%%% Calculate center of mass of array
    M = sum(expI);
    com = 0;
    for i=1:length(expI)
        com = expI(i)*i + com;
    end
    com;
    CoM = round(com/M);
end