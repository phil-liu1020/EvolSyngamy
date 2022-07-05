function mres0 = Invasion_Dynamics_Mass_Numerical_C_alphaprime(beta,C,alpha,A,M,tf,f0mut,theta1,deltam,deltaalpha,mres0,plots)

% This function iterates the invasion dynamics until fixation of the invading gamete is reached. In our model, we consider fixation to be reached at the time when the change in frequency of the invader from one generation to the next is less than theta1.
% List of parameters: deltam (and deltaalpha) - difference between mutant mass (mutant encounter rate) and resident gamete mass (mutant encounter rate).
%                     theta1 - change in frequency between a generation for the system to be deemed for fixation. This is usually a small parameter.
%                     tf - fertilisation period
%                     mres0 - resident gamete mass
%                     f0mut - initial frequency of mutant
%                     M - mass of adult
%                     A - number of adults at start of each generation
%                     alpha - encounter rate between gametes
%                     beta - resistance to survival of a gamete/agamete
%                     C - cost of fertilisation
%                     alpha is the encounter rate of the gamete (fixed)
%                     set plots to 1 if you want to plot the change in frequency of the mutant over the course of the invasion.
                
deltam=deltam*sign(randn);                                                       
g=0;                                                            
deltafmut=1;                                                               
i=0;

% g is the vector of mutant frequencies at the end of each fertilisation generation. It is initially set to 0 and grows in size.
% deltafmut is a parameter for the change in frequency of mutant at the end of each generation.
% Initially, deltafmut is set to a random large value (1 in this case) in order to start the while loop below.

% information for what happens inside the while loop:
%                     ode45 is applied to calculate the number of mutant
%                     and resident gametes at the end of the generation.                 
%                     exp(-beta/m) is S(beta,m), the survival probability of a body (either a resident or a mutant) of mass m.
%                     Parameters inside the loop:
%                                          Nres0 and Nmut0 are the numbers of resident and mutants in the gamete pool in the beginning of each fertilisation generation.
%                                          wmut and wres are the fitnesses of the mutant and resident respectively.
%                                          fmut is the frequency of the mutant at the each generation.
%                                          deltafmut is the change in mutant frequency at each generation. 

while deltafmut>theta1

   Nres0=(A*M*(1-f0mut))/mres0; 
   Nmut0=(A*M*f0mut)/(mres0+deltam);               
  
   [t,x] = ode45(@(t,x)unisexual_alphaprime_ODE(t,x,deltaalpha,alpha),[0,tf],[Nres0,Nmut0,0,0,0]);

   wmut = x(length(x(:,1)),5)*exp(-beta/(2*(mres0+deltam)))*(1-C) + x(length(x(:,1)),4)*exp(-beta/(2*mres0+deltam))*(1-C) + x(length(x(:,1)),2)*exp(-beta/(mres0+deltam)); 

   wres = x(length(x(:,1)),3)*exp(-beta/(2*mres0))*(1-C) + x(length(x(:,1)),4)*exp(-beta/(2*mres0+deltam))*(1-C) + x(length(x(:,1)),1)*exp(-beta/mres0);

   fmut=wmut/(wmut+wres);
   deltafmut=abs(fmut-f0mut);   
   
  
   f0mut=fmut;
   i = i + 1;

   if length(g)>10000
   break
   else
   end

   if i<=length(g)
   g(i)=f0mut;
   else
   g(end+1) = f0mut; 
   end

   
end
assignin('base','fmut',fmut)
assignin('base','i',i)
assignin('base','g',g)
assignin('base','deltam',deltam)

if plots==1
    plot(g)
    xlabel('t_g')
    ylabel('f_{mut}')
else
end

if fmut>0.5
    mres0 = mres0+deltam;
else
end

%%%%NOTE! as mutant macrogamete b can only invade if d/dfb(dfb/dt)>0
