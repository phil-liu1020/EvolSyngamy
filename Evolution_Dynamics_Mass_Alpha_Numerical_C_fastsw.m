function Evolution_Dynamics_Mass_Alpha_Numerical_C_fastsw(betaB,betaG,C,tsG,tsB,alpha0,A,M,tf,f0mut,fs0,theta1m,theta1a,deltam,deltaalpha,mres0,theta2,theta3)

% This function numerically plots the phase trajectory for m against \alpha, where m is the optimal size and \alpha the encounter rate that a
% gamete can evolve to, subject to environmental switching every fixed number of generations in each state. 
% tsG is the number of generations spent in the good state and tsB the number of generations in the bad state.

% List of parameters:
%                   tf - fertilisation period
%                   theta3 - small parameter for the change in mass for fixation.
%                   f0mut - initial frequnecy of gamete with mass m+\deltam, fs0 is the initial frequency of gamete with eocnunter rate \alpha'=\alpha+\delta\alpha.
%                   fs0 - initial frequency of gamete with encounter rate \alpha'=\alpha+\delta\alpha.
%                   theta1m is the theta1 for the gamete with mass m+\deltam and theta1a is the theta1 for the gamete with encounter rate alpha'
%                   mres0 is the initial mass of gamete
%                   If the trajectory changes by as small as theta3 in theta2 evolutionary time units \tau, then we assume the trajectrory has equilibrated.
%                   betaB and betaG are the resistance to survival in the bad and good environments respectiely.
%                   deltaalpha and deltam - mutational stepsize in \alpha and m respectively.
%
% Inside the while loop:
%                     mres and alphares - mass and encounter rate of resident
%                     B and G are binary variables representing the bad and good states respectively.
%                     g is the vector of mass, h the vector of encounter rate.
%                     If deltamres+deltaalphares changes by less than theta3 in theta2 invasion generations, then the
%                     system is deemed to have equilibrated.
%                     if sign(randn)=1, mutation occurs in the mass, else the muitation occurs in the encounter rate.
% Things to note:
% By default, we start in the good state.
% Mutations in mass and encounter rate do not occur simultaneously.

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



while abs(deltamres+deltaalphares)>theta3
tic

  
      if sign(randn)==1         
      mres=Invasion_Dynamics_Mass_Numerical_C_FastSwitchingEnv(betaB,betaG,tsG,tsB,C,alpha0,A,M,tf,f0mut,theta1m,deltam,mres0,0,B,G,0);  % Mutation in mass occurs, and if mutant invades, mres is the mutant mass, if mutant doesn't invade, mres is the resident mass.
      alphares=alpha0;  
   
      mres0=mres;
      g(end+1)=mres0;
      h(end+1)=alpha0;
     
      
      else
      alphares=Invasion_Dynamics_Alpha_Numerical_C_FastSwitchingEnv(betaB,betaG,tsG,tsB,C,alpha0,A,M,tf,fs0,theta1a,deltaalpha,mres0,0,B,G,0);    % if sign(randn)=-1, mutation occurs in \alpha the encounter rate
      mres=mres0;  
        
            if alpha0==deltaalpha
               break
            else
            end

      alpha0=alphares;  
      h(end+1)=alpha0;
      g(end+1)=mres0;
      end
                                                                                                       
  NINVGENS=NINVGENS+1;

  if mod(NINVGENS,theta2)==0                   % Checks whether the system has quasi equilibrated after theta2 invasions.
deltamres=abs(g(NINVGENS)-g(NINVGENS+1-theta2));
deltaalphares=abs(g(NINVGENS)-g(NINVGENS+1-theta2));
  else
  end



end




assignin('base','g',g)
assignin('base','h',h)