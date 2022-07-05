function FRTE_many_realisations_selectionstr(betaB,betaG,C,tsG,tsB,alpha0,A,M,tf,f0mut,fs0,theta1m,theta1a,deltam,deltaalpha,mres0,theta2,theta3,maxINVGENS,Nrealz)

% This code produces multiple realisations of stochastic trajectories under
% the FRTE regime i.e. multiple realisations of
% "Evolution_Dynamics_Mass_Alpha_Numerical_C_selectionstr_FRTE". This is done using a parfor loop.
% Nrealz is the number of stochastic trajectories we wish to obtain.


g=zeros(Nrealz,0);
h=zeros(Nrealz,0);


parfor i=1:Nrealz
[g(i,:),h(i,:)]=Evolution_Dynamics_Mass_Alpha_Numerical_C_selectionstr_slowsw3A(betaB,betaG,C,tsG,tsB,alpha0,A,M,tf,f0mut,fs0,theta1m,theta1a,deltam,deltaalpha,mres0,theta2,theta3,maxINVGENS);
end

assignin('base','g',g)
assignin('base','h',h)

save('gNrealz.mat','g')
save('hNrealz.mat','h')