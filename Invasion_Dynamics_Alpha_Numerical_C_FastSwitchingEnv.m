function alpha0 = Invasion_Dynamics_Alpha_Numerical_C_FastSwitchingEnv(betaB,betaG,tsG,tsB,C,alpha0,A,M,tf,f0mut,theta1,deltaalpha,m,plots,B,G, SCTE)


% This function simulates the invasion dynamics for a system with a
% mutation in \alpha, subject to a fast switching environment relative to the invasion timescale.

% List of parameters:
%                     deltaalpha - difference between mutant and resident encounter rate
%                     theta1 - change in mutant frequency between a fertilisation generation to be deemed for fixation.              
%                     tf - fertilisation period
%                     f0mut - initial frequency if mutant
%                     M - mass of adult
%                     A - number of adults at start of each generation
%                     alpha0 - initial encounter rate between gametes
%                     beta - resistance to survival of a gamete/agamete.
%                     C - cost of fertilisation
%                     B and G are the good and bad environmental conditions respectively.
%                     tsG and tsB are the times spent in the good and bad environments respectively.
%                     m is the mass of the gamete (fixed)
%                     Set SCTE to 1 if you want to assign a random sign to \deltaalpha.


% Information about the code:
% Nsw is the number of switches, defined by the number of jumps between the two for loops.
% tg represents which generation the system is in.
% The invasion dynamics are ran for a maximum of 10000 generations.

Nsw=1; 
tg=1; 
deltafmut=1;
if SCTE==1
else
deltaalpha=deltaalpha*sign(randn);  
end

g=zeros(1,50); 

                                                          
while deltafmut>theta1

    if B==1
    for i=1:tsB
    Nres0=(A*M*(1-f0mut))/m; 
    Nmut0=(A*M*f0mut)/m;

    [t,x] = ode45(@(t,x)unisexual_alphaprime_ODE(t,x,deltaalpha,alpha0),[0,tf],[Nres0,Nmut0,0,0,0]);
    wmut=x(length(x(:,1)),5)*exp(-betaB/(2*m))*(1-C)+x(length(x(:,1)),4)*exp(-betaB/(2*m))*(1-C)+x(length(x(:,1)),2)*exp(-betaB/m); 
    wres=x(length(x(:,1)),3)*exp(-betaB/(2*m))*(1-C)+x(length(x(:,1)),4)*exp(-betaB/(2*m))*(1-C)+x(length(x(:,1)),1)*exp(-betaB/m);
    fmut=wmut./(wmut+wres);
    deltafmut=abs(fmut-f0mut);
    if deltafmut<=theta1
    break
    else
    end
    f0mut=fmut;
    tg = tg + 1;
    if tg<=length(g)
    g(tg)=f0mut;
    else
    g(end+1) = f0mut; 
    end

    if length(g)>10000
    break
    else
    end

    end
    G=1;B=0;
    Nsw=Nsw+1;

    if length(g)>10000
    break
    else
    end

    elseif G==1
    
    for i=1:tsG
    Nres0=(A*M*(1-f0mut))/m; 
    Nmut0=(A*M*f0mut)/m;

    [t,x] = ode45(@(t,x)unisexual_alphaprime_ODE(t,x,deltaalpha,alpha0),[0,tf],[Nres0,Nmut0,0,0,0]);
    wmut=x(length(x(:,1)),5)*exp(-betaG/(2*m))*(1-C)+x(length(x(:,1)),4)*exp(-betaG/(2*m))*(1-C)+x(length(x(:,1)),2)*exp(-betaG/m); 
    wres=x(length(x(:,1)),3)*exp(-betaG/(2*m))*(1-C)+x(length(x(:,1)),4)*exp(-betaG/(2*m))*(1-C)+x(length(x(:,1)),1)*exp(-betaG/m);
    fmut=wmut./(wmut+wres);
    deltafmut=abs(fmut-f0mut); 
    if deltafmut<=theta1
    break
    else
    end
    f0mut=fmut;
    tg = tg + 1;
    if tg<=length(g)
    g(tg)=f0mut;
    else
    g(end+1) = f0mut; 
    end

    end
    Nsw=Nsw+1;

    if length(g)>10000
    break
    else
    end




    B=1; G=0;
    end

     
   
end
g(1)=[];
assignin('base','fmut',fmut)
assignin('base','g',g)
assignin('caller','deltaalpha',deltaalpha)
if Nsw==1
else
if B==1
    B=0; G=1;
else
    B=1; G=0;
end
end
assignin('caller','B',B)
assignin('caller','G',G)

if plots==1
    plot(g)
    
    xlabel('t_g')
    ylabel('f_{mut}')
else

end
assignin('caller','lengthg',length(g))

[Max,Imax] = max(g);
[Min,Imin] = min(g);
if Imax-Imin>0
    alpha0 = alpha0+deltaalpha;
else
end
end