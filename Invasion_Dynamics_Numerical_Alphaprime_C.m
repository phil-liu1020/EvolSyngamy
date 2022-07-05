function alphares0 =Invasion_Dynamics_Numerical_Alphaprime_C(beta,C,alphares0,A,M,tf,f0mut,theta1,deltaalpha,m,plots)

% This function iterates the invasion dynamics until fixation of the invading gamete with different encounter rate is reached. In our model, we consider fixation to be reached at the time when the change in frequency of the invader from one generation to the next is less than theta1.
% List of parameters: deltaalpha - difference between mutant encounter rate and resident encounter rate.
%                     theta1 - change in frequency between a generation for the system to be deemed for fixation. This is usually a small parameter.
%                     tf - fertilisation period
%                     alphares0 - resident gamete encounter rate
%                     f0mut - initial frequency of mutant
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
%                                          Nres0 and Nmut0 are the numbers of resident and mutants in the gamete pool in the beginning of each fertilisation generation.
%                                          ws and wa are the fitnesses of the mutant and resident respectively.
%                                          fmut is the frequency of the mutant at the each generation.
%                                          deltaalpha is the change in mutant frequency at each generation. 



if alphares0-deltaalpha<=0                                                     
else
deltaalpha=deltaalpha*sign(randn);                                         
end
g=0;                                                             
deltafmut=1;                                                                 
i=0;

while deltafmut>theta1

   Nres0=(A*M*(1-f0mut))/m;                                                    
   Nmut0=(A*M*f0mut)/m;                                                        

   [t,x] = ode45(@(t,x)unisexual_alphaprime_ODE(t,x,deltaalpha,alphares0),[0,tf],[Nres0,Nmut0,0,0,0]);  

   ws=(Nmut0-x(length(x(:,1)),2))*exp(-beta/(2*m))*(1-C)+x(length(x(:,1)),2)*exp(-beta/m);    

   wa=(Nres0-x(length(x(:,1)),1))*exp(-beta/(2*m))*(1-C)+x(length(x(:,1)),1)*exp(-beta/m);

   fmut=ws/(ws+wa);
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

if plots==1
    plot(g)
    xlabel('t_g')
    ylabel('f_{mut}')
else
end

if fmut>0.5
    alphares0 = alphares0+deltaalpha;
else
end

assignin('base','deltaalpha',deltaalpha)