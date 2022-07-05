function Evolution_Dynamics_Mass_Alpha_Numerical_C(beta,C,alpha0,A,M,tf,f0mut,fs0,theta1m,theta1a,deltam,deltaalpha,mres0,theta2,theta3)

% This function produces a stochastic trajectory to show the mass and \alpha that a gamete can coevolve to, given particular initial conditions on \alpha and m.
% If the trajectory changes by as small as theta3 in theta2 evolutionary time units \tau, then we assume the trajectory has equilibrated.
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

g=zeros(1,theta2);  
h=zeros(1,theta2); 
deltamres=1;
deltaalphares=1;
NLOOPS=0;


% g and h are vectors for the mass and \alpha respectively (use the command "plot(g,h)" to see the numerical trajectory on the \alpha-m plane). Their lengths are initially set to theta2, and grows as more generations are elapsed.
% deltafmut and deltaalphares are both set to a random large value of 1 in order to start the while loop below.

% Notes for things happening inside the while loop:
%                                              Mutations in mass and encounter rate do not occur simultaneously.
%                                              If sign(randn)=1, mutation occurs in the mass, otherwise mutations occur in \alpha.
%                                              deltaalpha is set to 0 in the function 
%                                              Parameters inside the loop:
%                                                                  mres and alphares - mass and encounter rate of the gamete at each invasion generation respectively
%                                                                  Nmut - the number of random mutation events that've occured
%                                                                  deltamres and deltaalphares - change in mass and encounter rate respectively after theta2 generations
%                                                                  NLOOPS - number of theta2 generations that've elapsed.
                                                           

while abs(deltamres+deltaalphares)>theta3
tic
  for Nmut=1:theta2
      if sign(randn)==1                                                                                          
      mres=Invasion_Dynamics_Mass_Numerical_C_alphaprime(beta,C,alpha0,A,M,tf,f0mut,theta1m,deltam,0,mres0,0);  
      alphares=alpha0;                                                                                           
      else
      alphares=Invasion_Dynamics_Numerical_Alphaprime_C(beta,C,alpha0,A,M,tf,fs0,theta1a,deltaalpha,mres0,0);    
      mres=mres0;
      end                                                                                                        

  alpha0=alphares;
  mres0=mres;

  g(Nmut+theta2*NLOOPS)=mres0;
  h(Nmut+theta2*NLOOPS)=alpha0;
  end

deltamres=g(Nmut+theta2*NLOOPS)-g(theta2*NLOOPS+1);
deltaalphares=g(Nmut+theta2*NLOOPS)-g(theta2*NLOOPS+1);

NLOOPS=NLOOPS+1;

fprintf('Processing %d...',NLOOPS);


if alpha0==deltaalpha
    break
else
end

toc
end

assignin('base','g',g)
assignin('base','h',h)