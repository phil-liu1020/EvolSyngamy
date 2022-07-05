function Evolution_Dynamics_Mass_Numerical_C(beta,C,alpha,A,M,tf,f0mut,theta1,deltam,mres0,asexualODE,theta2,theta3,plots)

% This function numerically evaluates the optimal size that a gamete can evolve to, given particular initial conditions on \alpha and m.
% If the mass changes by as small as theta3 in theta2 evolutionary time units \tau, then we assume the mass has equilibrated.
% List of parameters:
%                     deltam - difference between mutant and resident mass
%                     theta1 - change in frequency between a fertilisation generation to be deemed for fixation. 
%                     theta2 - number of invasion generations after which the code checks the change in mass. This is usually of the order 10.
%                     theta3 - change in mass between theta2 invasion generations to be deemed for fixation. 
%                     tf - fertilisation period
%                     mres0 - resident gamete mass
%                     f0mut - initial freqneucy if mutant
%                     M - mass of adult
%                     A - number of adults at start of each generation
%                     alpha - encounter rate between gametes
%                     beta - resistance to survival of a gamete/agamete.
%                     C - cost of fertilisation
%                     Set asexualODE to 1 to investigate the system 'asexual_ODE', and to 0 to investigate the system 'unisexual_ODE' where mutant unisexual invades resident unisexual.

g=zeros(1,theta2);  
deltamres=1;
NLOOPS=0;

% g is a vector of mutant frequencies at the end of each fertilisation generation. Its length is initially set to [0,theta2] and grows in size as more generations are elapsed.(use the command "plot(g)" to see how the mass evolves over time)
% deltamres is set to a random large value of 1 in order to start the while loop below.
% If |deltamres|<theta3, then the system is deemed to have equilibrated
% deltamres is set to a random large value of 1 in order to start the while loop below.
% Parameters inside the loop:
%                                                                  mres- mass of the gamete at each invasion generation
%                                                                  Nmut - the number of random mutation events that've occured
%                                                                  deltamres - change in mass after theta2 generations
%                                                                  NLOOPS - number of theta2 generations that've elapsed.


while abs(deltamres)>theta3
tic
  for Nmut=1:theta2 
  mres=Invasion_Dynamics_Mass_Numerical_C(beta,C,alpha,A,M,tf,f0mut,theta1,deltam,mres0,0,asexualODE);
  mres0=mres;
  g(Nmut+theta2*NLOOPS)=mres0;
  end

deltamres=g(Nmut+theta2*NLOOPS)-g(theta2*NLOOPS+1);
NLOOPS=NLOOPS+1;

fprintf('Processing %d...',NLOOPS);
toc
end

assignin('base','g',g)

if plots==1
    plot(g)
    xlabel('\tau')
    ylabel('m')
else
end