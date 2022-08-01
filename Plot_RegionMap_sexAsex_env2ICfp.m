function Plot_RegionMap_sexAsex_env2ICfp(C,A,M,T,k)

% This checks what \alpha the system tends to if starting from the
% switching induced fixed point. If the switching induced fixed point
% doesn't exist, beta2AS(j,i)=-1.

% beta2AS(j,i)=0 if asexuality selected for
% beta2AS(j,i)=1 if sexuality selected for
% k is the proportion of time the system spends in the good state. Also defined as lambdaBG/(lambdaBG+lambdaGB)



b1=(0.005:0.005:4);
b2=(0.005:0.005:4);    




[beta1,beta2] = meshgrid(b1,b2);
beta2AS=zeros(length(b2),length(b1)); 

alphastarNUM=(1 + (C-1)*exp (2*beta1./(beta2.*(1 - k) + beta1.*k)).*k + (C - 1).*(1 - k).*exp (2*beta2./(beta2.*(1 - k) + beta1.*k))).*(beta1.*k + beta2.*(1 - k));
alphastarDEN=4*A*(C-1)*M*T*(k.*exp (2*beta2./(beta2.*(1 - k) + beta1.*k)) + exp(2*beta1./(beta2.*(1 - k) + beta1.*k)).*((1 - k) + (C - 1)*exp (2*beta2./(beta2.*(1 - k) + beta1.*k))) );
alphastar=alphastarNUM./alphastarDEN;
alphastar=alphastar';
mstar=0.25*(max(beta1,beta2).*(1-k) + min(beta1,beta2).*k);

   
for i=1:length(b1)
    for j=1:length(b2)
        if alphastar(j,i)>0
        phasetraj_alphaprime(A,M,T,b2(i),C,alphastar(j,i),mstar(j,i),0.1,0)
           if alphaend<0.000001
           beta2AS(j,i)=0;    
           else
           beta2AS(j,i)=1;      
           end
        else
        beta2AS(j,i)=-1;    
        end
    
      
    end
end


assignin('base','beta2AS',beta2AS)
assignin('base','mstar',mstar)
assignin('base','alphastar',alphastar)