function dx = unisexual_alphaprime_ODE(t,x,deltaalpha,alpha)

% This is the ODE for the fertilisation kinetics in the presence of mutant with different encounter rate \alpha. The gametes are assumed
% to be unisexual meaning that all gametes have the ability to fertilise with one another. NB this function can also be used to describe the fertilisation kinetics of a mutant with different mass, but in this case set \deltaalpha=0. 

%x(1) is the number of resident gametes
%x(2) is the number of mutant gametes
%x(3) is the number of zygotes that are "resident-resident complexes"
%x(4) is the number of zygotes that are "resident-mutant complexes"
%x(5) is the number of zygotes that are "mutant-mutant complexes"

% The resident fertilises with  mutant at encounter rate (\alpha+\alpha')/2
% where \alpha'=\alpha+\deltaalpha.

dx = zeros(2,1);

dx(1) = -alpha*x(1)^2-0.5*(2*alpha+deltaalpha)*x(2)*x(1);

dx(2) = -(alpha+deltaalpha)*x(2)^2-0.5*(2*alpha+deltaalpha)*x(2)*x(1);

dx(3) = alpha*x(1)^2;

dx(4) = 0.5*(2*alpha+deltaalpha)*x(1)*x(2);   

dx(5) = (alpha+deltaalpha)*x(2)^2;  
