function Plot_RegionMap_sexAsex_env1ICfp(C,A,M,T,betarange,k)

% This checks what \alpha the system tends to if starting from the
% switching induced fixed point. If the switching induced fixed point
% doesn't exist, beta1AS(j,i)=-1.

% beta1AS(j,i)=0 if asexuality selected for
% beta1AS(j,i)=1 if sexuality selected for
% k is the proportion of time the system spends in the good state. Also defined as lambdaBG/(lambdaBG+lambdaGB)
% Set betarange to 1 if plotting the zoomed in version
% Set betarange to 0 if plotting zoomed out version

if betarange==1
b1=(0.002:0.002:0.2);
b2=(3.01:0.01:4);
else
b1=(0.005:0.005:4);
b2=(0.005:0.005:4);    
end


[beta1,beta2] = meshgrid(b1,b2);
beta1AS=zeros(length(b2),length(b1)); 

alphastarNUM=(1 + (C-1)*exp (2*beta1./(beta2.*(1 - k) + beta1.*k)).*k + (C - 1).*(1 - k).*exp (2*beta2./(beta2.*(1 - k) + beta1.*k))).*(beta1.*k + beta2.*(1 - k));
alphastarDEN=4*A*(C-1)*M*T*(k.*exp (2*beta2./(beta2.*(1 - k) + beta1.*k)) + exp(2*beta1./(beta2.*(1 - k) + beta1.*k)).*((1 - k) + (C - 1)*exp (2*beta2./(beta2.*(1 - k) + beta1.*k))) );
alphastar=alphastarNUM./alphastarDEN;
alphastar=alphastar';
mstar=0.25*(max(beta1,beta2).*(1-k) + min(beta1,beta2).*k);
mstar=mstar';
   
for i=1:length(b1)
    for j=1:length(b2)
        if alphastar(j,i)>0
        phasetraj_alphaprime(A,M,T,b1(j),C,alphastar(j,i),mstar(j,i),0.1,0)
           if alphaend<0.000001
           beta1AS(j,i)=0;    
           elseif alphaend-alphastar(j,i)>0
           beta1AS(j,i)=1;
           else
           beta1AS(j,i)=10;   
           end
        else
        beta1AS(j,i)=-1;    
        end
    
      
    end
end


assignin('base','beta1AS',beta1AS)
assignin('base','mstar',mstar)
assignin('base','alphastar',alphastar)

% imagesc(0.002:0.002:0.2,3.01:0.01:4,sumAS')
% set(gca,'YDir','normal')