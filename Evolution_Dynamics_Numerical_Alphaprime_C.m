function Evolution_Dynamics_Numerical_Alphaprime_C(beta,C,alphares0,A,M,tf,f0mut,theta1,deltaalpha,m,theta2,theta3,plots)

% This function numerically evaluates the optimal size that a gamete can evolve to, given particular initial conditions on \alpha and m.
% If the trajectory changes by as small as theta3 in theta2 evolutionary time units \tau, then we assume the trajectory has equilibrated.
% List of parameters:
%                     deltam - difference between mutant and resident mass
%                     theta1 - change in frequency between a fertilisation generation to be deemed for fixation. 
%                     theta2 - number of invasion generations after which the code checks the change in encounter rate. theta2 is usually of the order 10.
%                     theta3 - change in mass between theta2 invasion generations to be deemed for fixation. 
%                     tf - fertilisation period
%                     mres0 - resident gamete mass
%                     f0mut - initial freqneucy if mutant
%                     M - mass of adult
%                     A - number of adults at start of each generation
%                     alphares0 - initial encounter rate between gametes
%                     beta - resistance to survival of a gamete/agamete.
%                     C - cost of fertilisation
%                     Set asexualODE to 1 to investigate the system 'asexual_ODE', and to 0 to investigate the system 'unisexual_ODE' where mutant unisexual invades resident unisexual.
% set plots=1 if you want to see the phase plane plotted out.

g=zeros(1,theta2);  
deltaalphares=1;
NLOOPS=0;

% g is a vector of mutant frequencies at the end of each fertilisation generation. Its length is initially set to theta2 and grows in size as more generations are elapsed.(use the command "plot(g)" to see how \alpha evolves over time)
% deltaalphares is set to a random large value of 1 in order to start the while loop below.
% If |deltaalphares|<theta3, then the system is deemed to have equilibrated.
% deltaalphares is set to a random large value of 1 in order to start the while loop below.
% Parameters inside the loop:
%                                                                  alphares- encounter rate of the gamete at each invasion generation
%                                                                  Nmut - the number of random mutation events that've occured
%                                                                  deltaalphares - change in mass after theta2 generations
%                                                                  NLOOPS - number of theta2 generations that've elapsed.

while abs(deltaalphares)>theta3
tic
  for Nmut=1:theta2
  alphares=Invasion_Dynamics_Numerical_Alphaprime_C(beta,C,alphares0,A,M,tf,f0mut,theta1,deltaalpha,m,0);
  alphares0=alphares;
  g(Nmut+theta2*NLOOPS)=alphares0;
  end

deltaalphares=g(Nmut+theta2*NLOOPS)-g(theta2*NLOOPS+1);
NLOOPS=NLOOPS+1;

fprintf('Processing %d...',NLOOPS);
toc
end

assignin('base','g',g)

if plots==1
    plot(g)
    xlabel('\tau')
    ylabel('\alpha')
else
end