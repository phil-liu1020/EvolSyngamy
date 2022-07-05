function Evolution_Dynamics_Mass_Alpha_Numerical_C_fixationprob_fastsw(betaB,betaG,C,tsG,tsB,alphares0,A,M,tf,f0mut,fs0,theta1m,theta1a,deltam,deltaalpha,mres0,theta2,theta3)

% This code simulates the coevolutionary dynamics for \alpha and m in an environment that switches after a fixed time period. The switching rate is fast relative to the invasion timescale (FRTI). No plastic mechanism is assumed here, so the dynamics undergo bet-hedging.
% If the trajectory changes by as small as theta3 in theta2 units of evolutionary time \tau, then we assume the trajectory has equilibrated. alphares0 and mres0 are initial values of the encounter rate and mass,
% where we want to start our phase trajectory.

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
%                     alphares0 - initial encounter rate between gametes
%                     beta - resistance to survival of a gamete/agamete.
%                     C - cost of fertilisation
%                     B and G represent the states this system are in, B for bad state and G for good state.
%                     tsG and tsB are the amount of times spent in the good and bad states respectively.

% Parameters inside the while loop:
%                                                                  alphares- encounter rate of the gamete at each invasion generation
%                                                                  deltaalphares - change in mass after theta2 generations
%                                                                  NINVGENS - number of invasion generations elapsed.


g=[];  
h=[]; 
deltamres=1;
deltaalphares=1;

if tsG==0
G=0; B=1;                        
elseif tsB==0
G=1; B=0;
else
G=1; B=0;                         
end

NINVGENS=0;

% If no time is spent in the good state, then the system starts in the bad state. By default, we start in the good state.
% g and h are vectors of mutant frequencies at the end of each
% fertilisation generation. Its length is initially set to 0 and grows in size as more generations are elapsed.(use the command "plot(g)" to see how \alpha evolves over time).
% deltaalphares  and deltamres are set to a random large value of 1 in order to start the while loop below.
% If |deltaalphares|<theta3, then the system is deemed to have equilibrated.
% Parameters inside the loop:
%                          alphares- encounter rate of the gamete at each invasion generation
%



while abs(deltamres+deltaalphares)>theta3
tic

  
      if sign(randn)==1         % if sign(randn)=1, mutation occurs in the mass
      mres=Invasion_Dynamics_Mass_Numerical_C_FastSwitchingEnv(betaB,betaG,tsG,tsB,C,alphares0,A,M,tf,f0mut,theta1m,deltaalpha,mres0,0,B,G,0);  % Mutation in mass occurs, and if mutant invades, mres is the mutant mass, if mutant doesn't invade, mres is the resident mass.
      alphares=alphares0;                                                                                           % \alpha doesn't change here, hence alphaprime is set to 0 in the function in the line above.
    
      if mres~=mres0            % if the mutation in mass has occured, then we let fixation occur with probability pfix. 
        s=abs(((1/tsB)/(1/tsB+1/tsG))*(-(4*mres0.*(mres0-betaG)-A*(C-1).*exp(betaG./(2*mres0))*M.*(4*mres0-betaG).*alphares0*tf)./(4*mres0.^2.*(mres0-A*M.*alphares0*tf*(C-1).*exp(betaG./(2*mres0))))*deltam) + ((1/tsG)/(1/tsB+1/tsG))*(-(4*mres0.*(mres0-betaB)-A*(C-1).*exp(betaB./(2*mres0))*M.*(4*mres0-betaB).*alphares0*tf)./(4*mres0.^2.*(mres0-A*M.*alphares0*tf*(C-1).*exp(betaB./(2*mres0))))*deltam));
        N=(A*M*(1-f0mut))/mres0+(A*M*f0mut)/(mres0+deltam); 
        pfix=(1-exp(-4*s*N*f0mut))/(1-exp(-4*s*N));
          if rand<pfix                              
             mres0=mres;        % the new resident mass is the mutant mass if fixation has occured
             g(end+1)=mres0;
             h(end+1)=alphares0;
          else
             g(end+1)=mres0;
             h(end+1)=alphares0;
          end
  
     else
 
     mres0=mres;               % the mutant is of the same mass as the resident, hence this line deoesn't do anything.
     g(end+1)=mres0;
     h(end+1)=alphares0;
     end
      
      else
      alphares=Invasion_Dynamics_Alpha_Numerical_C_FastSwitchingEnv(betaB,betaG,tsG,tsB,C,alphares0,A,M,tf,fs0,theta1a,deltaalpha,mres0,0,B,G,0);    % if sign(randn)=-1, mutation occurs in \alpha the encounter rate
      mres=mres0;  
        if alphares~=alphares0      % if the mutation in alpha has occured, then we let fixation occur with probability pfix. 
           s=abs(((1/tsB)/(1/tsB+1/tsG))*( ((1+(C-1)*exp(betaG./(2*mres0))).*mres0.*log(1+(A*M.*alphares0*tf)./mres0))./(-2.*mres0.*alphares0+2*A*(C-1).*exp(betaG./(2*mres0))*M.*alphares0*tf.*alphares0) )*deltaalpha + ((1/tsG)/(1/tsB+1/tsG))* ( ((1+(C-1)*exp(betaB./(2*mres0))).*mres0.*log(1+(A*M.*alphares0*tf)./mres0))./(-2.*mres0.*alphares0+2*A*(C-1).*exp(betaB./(2*mres0))*M.*alphares0*tf.*alphares0) )*deltaalpha);
           N=(A*M*(1-fs0))/mres0+(A*M*fs0)/mres0; 
           pfix=(1-exp(-4*s*N*fs0))/(1-exp(-4*s*N));
            if rand<pfix
               alphares0=alphares;     % the new resident alpha is the mutant alpha if fixation has occured
               h(end+1)=alphares0;
               g(end+1)=mres0;
            else
               h(end+1)=alphares0;
               g(end+1)=mres0;
            end

            if alphares0==deltaalpha
               break
            else
            end
  
       else
 
       alphares0=alphares;            % the mutant is of the same alpha as the resident, hence this line doesn't do anything
       h(end+1)=alphares0;
       g(end+1)=mres0;
       end
      end                                                                                                        % Mutations in mass and encounter rate do not occur simultaneously .
  NINVGENS=NINVGENS+1;

  if mod(NINVGENS,theta2)==0                   % Checks whether the system has quasi equilibrated after theta2 invasions.
deltamres=abs(g(NINVGENS)-g(NINVGENS+1-theta2));
deltaalphares=abs(g(NINVGENS)-g(NINVGENS+1-theta2));
  else
  end



end




assignin('base','g',g)
assignin('base','h',h)