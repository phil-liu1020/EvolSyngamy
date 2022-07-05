function alpha0 =Invasion_Dynamics_Numerical_Alphaprime_C(beta,C,alpha0,A,M,tf,fs0,theta1,deltaalpha,m,plots)

% This function iterates the invasion dynamics until fixation of the invading gamete with different encounter rate is reached. In our model, we consider fixation to be reached at the time when the change in frequency of the invader from one generation to the next is less than theta1.
% List of parameters: deltaalpha - difference between mutant encounter rate and resident encounter rate.
%                     theta1 - change in frequency between a generation for the system to be deemed for fixation. This is usually a small parameter.
%                     tf - fertilisation period
%                     alpha0 - resident gamete encounter rate
%                     fs0 - initial frequency of mutant
%                     M - mass of adult
%                     A - number of adults at start of each generation
%                     beta - resistance to survival of a gamete/agamete
%                     C - cost of fertilisation
%                     m - gamete mass (fixed)
%                     set plots to 1 if you want to plot the change in frequency of the mutant over the course of the invasion.


% g is the vector of mutant frequencies at the end of each fertilisation generation. It is initially set to 0 and grows in size.
% deltafmut is a parameter for the change in frequency of mutant at the end of each generation.
% Initially, deltafmut is set to a random large value (1 in this case) in order to start the while loop below.

% information for what happens inside the while loop:
%                     ode45 is applied to calculate the number of mutant
%                     and resident gametes at the end of the generation.                 
%                     exp(-beta/m) is S(beta,m), the survival probability of a body (either a resident or a mutant) of mass m.
%                     if the change in encounter rate is negative and the initial encounter minus the change in encounter rate is less than 0, then deltaalpha is what it is, otherwise multiply deltaalpha by a randomly generated sign.
%                     Parameters inside the loop:
%                                          Na0 and Ns0 are the numbers of resident and mutants in the gamete pool in the beginning of each fertilisation generation.
%                                          ws and wa are the fitnesses of the mutant and resident respectively.
%                                          fs is the frequency of the mutant at the each generation.
%                                          deltaalpha is the change in mutant frequency at each generation. 



if alpha0-deltaalpha<=0                                                     
else
deltaalpha=deltaalpha*sign(randn);                                         
end
g=0;                                                             
deltafs=1;                                                                 
i=0;

while deltafs>theta1

   Na0=(A*M*(1-fs0))/m;                                                    
   Ns0=(A*M*fs0)/m;                                                        

   [t,x] = ode45(@(t,x)unisexual_alphaprime_ODE(t,x,deltaalpha,alpha0),[0,tf],[Na0,Ns0,0,0,0]);  

   ws=(Ns0-x(length(x(:,1)),2))*exp(-beta/(2*m))*(1-C)+x(length(x(:,1)),2)*exp(-beta/m);    

   wa=(Na0-x(length(x(:,1)),1))*exp(-beta/(2*m))*(1-C)+x(length(x(:,1)),1)*exp(-beta/m);

   fs=ws/(ws+wa);
   deltafs=abs(fs-fs0);                                          
   
  
   fs0=fs;
   i = i + 1;

   if length(g)>10000
   break
   else
   end

   if i<=length(g)
   g(i)=fs0;
   else
   g(end+1) = fs0; 
   end
   
end
assignin('base','fs',fs)
assignin('base','i',i)
assignin('base','g',g)

if plots==1
    plot(g)
    xlabel('t_g')
    ylabel('f_{mut}')
else
end

if fs>0.5
    alpha0 = alpha0+deltaalpha;
else
end

assignin('base','deltaalpha',deltaalpha)