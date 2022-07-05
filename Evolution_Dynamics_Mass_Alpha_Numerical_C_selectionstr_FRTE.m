function [g,h]=Evolution_Dynamics_Mass_Alpha_Numerical_C_selectionstr_FRTE(betaB,betaG,C,tsG,tsB,alpha0,A,M,tf,f0mut,fs0,theta1m,theta1a,deltam,deltaalpha,mres0,theta2,theta3,maxINVGENS)

% This code simulates the coevolutionary dynamics for \alpha and m in an
% environment that switches after a fixed time period. No plastic mechanism is assumed here, so the dynamics undergo bet-hedging.

% This code accounts for the selection strengths in a switching environment. By accounting for the selection
% strengths, the numerically simulated trajectories will follow the analytical trajectories closer in the phase plane.
% If the trajectory changes by as small as theta3 in theta2 units of evolutionary time \tau, then we assume the trajectory has equilibrated.
% alpha0 and mres0 are initial values of encounter rate and mass, where we want to start our phase trajectory



% List of parameters:
%                     deltam - difference between mutant and resident mass
%                     deltaalpha - difference between mutant and resident encounter rate
%                     theta1 - change in frequency between a fertilisation generation to be deemed for fixation. 
%                         (theta1m - for a mutant with different mass)
%                         (theta1a - for a mutant with different encounter rate)
%                     theta2 - number of invasion generations after which the code checks the change in mass. This is usually of the order 10.
%                     theta3 - change in mass between theta2 invasion generations to be deemed for fixation. 
%                     tf - fertilisation period
%                     mres0 - resident gamete mass
%                     f0mut - initial freqneucy if mutant with different mass
%                     fs0 - initial freqneucy if mutant with different encounter rate
%                     M - mass of adult
%                     A - number of adults at start of each generation
%                     alpha0 - initial encounter rate between gametes
%                     beta - resistance to survival of a gamete/agamete.
%                     C - cost of fertilisation
%                     B and G represent the states this system are in, B for bad state and G for good state.
%                     tsG and tsB are the amount of times spent in the good and bad states respectively.
%                     maxINVGENS - maximum number of generations to run this function for.

% Parameters inside the while loop:
%                                                                  alphares- encounter rate of the gamete at each invasion generation
%                                                                  deltaalphares - change in mass after theta2 generations
%                                                                  NINVGENS (g and b) - number of invasion generations elapsed (in the good and bad environments).
%                                                                  Sm and Salpha (G and B) - selection strength of a mutant (with a different mass and encounter rate respectively)
%                                                                  pfix - The probability that a given mutant will fixate.
%                                                                
%                                                                  
%                                                                

                           
g=[];  
h=[]; 
deltamres=1;
deltaalphares=1;

if tsG==0
G=0; B=1;                        % If no time is spent in the good state, then the system starts in the bad state
elseif tsB==0
G=1; B=0;
else
G=1; B=0;                         % By default, we start in the good state.
end

NINVGENS=0;
NINVGENSg=0;
NINVGENSb=0;


while abs(deltamres+deltaalphares)>theta3
tic
  

if G==1
NINVGENSb=0;

      smG=abs(-(4*mres0.*(mres0-betaG )-A*(C-1).*exp(betaG./(2*mres0))*M.*(4*mres0-betaG).*alpha0*tf)./(4*mres0.^2.*(mres0-A*M.*alpha0*tf*(C-1).*exp(betaG./(2*mres0))))*deltam);
      salphaG=abs( ( ((1+(C-1)*exp(betaG./(2*mres0))).*mres0.*log(1+(A*M.*alpha0*tf)./mres0))./(-2.*mres0.*alpha0+2*A*(C-1).*exp(betaG./(2*mres0))*M.*alpha0*tf.*alpha0) )*deltaalpha);
      smB=abs(-(4*mres0.*(mres0-betaB )-A*(C-1).*exp(betaB./(2*mres0))*M.*(4*mres0-betaB).*alpha0*tf)./(4*mres0.^2.*(mres0-A*M.*alpha0*tf*(C-1).*exp(betaB./(2*mres0))))*deltam);
      salphaB=abs( ( ((1+(C-1)*exp(betaB./(2*mres0))).*mres0.*log(1+(A*M.*alpha0*tf)./mres0))./(-2.*mres0.*alpha0+2*A*(C-1).*exp(betaB./(2*mres0))*M.*alpha0*tf.*alpha0) )*deltaalpha);


      if rand<smG/(smG+salphaG)   % mass mutation chosen with probability equal to the ratio of selection strength in m to selection strength in \alpha.                                                                                
         mres=Invasion_Dynamics_Mass_Numerical_C_alphaprime(betaG,C,alpha0,A,M,tf,f0mut,theta1m,deltam,0,mres0,0);
         alphares=alpha0; 

         if mres~= mres0
         pfix=smG/(smG+smB);
          if rand<pfix                              
             mres0=mres;        
             g(end+1)=mres0;
             h(end+1)=alpha0;
          else
             g(end+1)=mres0;
             h(end+1)=alpha0;
          end         
         else
         mres0=mres;              
         g(end+1)=mres0;
         h(end+1)=alpha0;

         end


      else    % alpha mutation
         alphares=Invasion_Dynamics_Numerical_Alphaprime_C(betaG,C,alpha0,A,M,tf,fs0,theta1a,deltaalpha,mres0,0);
         mres=mres0;

         if alphares~= alpha0
         pfix=salphaG/(salphaG+salphaB);
          if rand<pfix                              
             alpha0=alphares;        
             g(end+1)=mres0;
             h(end+1)=alpha0;
          else
             g(end+1)=mres0;
             h(end+1)=alpha0;
          end         
         else
         alpha0=alphares;              
         g(end+1)=mres0;
         h(end+1)=alpha0;

         end        
      end
 
  
NINVGENS=NINVGENS+1;
NINVGENSg=NINVGENSg+1;

if NINVGENS==maxINVGENS
break
else
end


  if mod(NINVGENSg,tsG)==0 && tsB~=0           % Checks whether we need to switch from good to bad state.
      B=1; G=0;
  else
  end

  if mod(NINVGENS,theta2)==0                   % Checks whether the system has quasi equilibrated after theta2 invasions.
deltamres=abs(g(NINVGENS)-g(NINVGENS+1-theta2));
deltaalphares=abs(h(NINVGENS)-h(NINVGENS+1-theta2));
  else
  end



elseif B==1
NINVGENSg=0;

      smG=abs(-(4*mres0.*(mres0-betaG )-A*(C-1).*exp(betaG./(2*mres0))*M.*(4*mres0-betaG).*alpha0*tf)./(4*mres0.^2.*(mres0-A*M.*alpha0*tf*(C-1).*exp(betaG./(2*mres0))))*deltam);
      salphaG=abs( ( ((1+(C-1)*exp(betaG./(2*mres0))).*mres0.*log(1+(A*M.*alpha0*tf)./mres0))./(-2.*mres0.*alpha0+2*A*(C-1).*exp(betaG./(2*mres0))*M.*alpha0*tf.*alpha0) )*deltaalpha);
      smB=abs(-(4*mres0.*(mres0-betaB )-A*(C-1).*exp(betaB./(2*mres0))*M.*(4*mres0-betaB).*alpha0*tf)./(4*mres0.^2.*(mres0-A*M.*alpha0*tf*(C-1).*exp(betaB./(2*mres0))))*deltam);
      salphaB=abs( ( ((1+(C-1)*exp(betaB./(2*mres0))).*mres0.*log(1+(A*M.*alpha0*tf)./mres0))./(-2.*mres0.*alpha0+2*A*(C-1).*exp(betaB./(2*mres0))*M.*alpha0*tf.*alpha0) )*deltaalpha);
      
      if rand<smB/(smB+salphaB)   % mass mutation                                                                                     
         mres=Invasion_Dynamics_Mass_Numerical_C_alphaprime(betaB,C,alpha0,A,M,tf,f0mut,theta1m,deltam,0,mres0,0);
         alphares=alpha0; 

         if mres~= mres0
         pfix=smB/(smG+smB);
          if rand<pfix                              
             mres0=mres;        
             g(end+1)=mres0;
             h(end+1)=alpha0;
          else
             g(end+1)=mres0;
             h(end+1)=alpha0;
          end         
         else
         mres0=mres;              
         g(end+1)=mres0;
         h(end+1)=alpha0;

         end


      else    % alpha mutation
         alphares=Invasion_Dynamics_Numerical_Alphaprime_C(betaB,C,alpha0,A,M,tf,fs0,theta1a,deltaalpha,mres0,0);
         mres=mres0;

         if alphares~= alpha0
         pfix=salphaB/(salphaG+salphaB);
          if rand<pfix                              
             alpha0=alphares;        
             g(end+1)=mres0;
             h(end+1)=alpha0;
          else
             g(end+1)=mres0;
             h(end+1)=alpha0;
          end         
         else
         alpha0=alphares;              
         g(end+1)=mres0;
         h(end+1)=alpha0;

         end        
      end

  NINVGENS=NINVGENS+1;
  NINVGENSb=NINVGENSb+1;

  if NINVGENS==maxINVGENS
  break
  else
  end

  
  if mod(NINVGENSb,tsB)==0 && tsG~=0
      B=0; G=1;
  else
  end

  if mod(NINVGENS,theta2)==0
deltamres=abs(g(NINVGENS)-g(NINVGENS+1-theta2));
deltaalphares=abs(h(NINVGENS)-h(NINVGENS+1-theta2));
  else
  end


end

end


assignin('base','g',g)
assignin('base','h',h)

end


