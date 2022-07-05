function dx = unisexual_ODE(t,x,alpha)

% This is the ODE for the fertilisation kinetics in the presence of mutant of a different mass. The gametes are assumed
% to be unisexual, meaning that all gametes have the ability to fertilise with one another.

%x(1) is the number of resident gametes
%x(2) is the number of mutant gametes
%x(3) is the number of zygotes that are "resident-resident complexes"
%x(4) is the number of zygotes that are "resident-mutant complexes"
%x(5) is the number of zygotes that are "mutant-mutant complexes"

% The mutant has mass m+\deltam.

dx = zeros(5,1);

dx(1) = -alpha*(x(1)^2+x(1)*x(2));

dx(2) = -alpha*(x(1)*x(2)+x(2)^2);

dx(3) = alpha*(x(1)^2);   

dx(4) = alpha*x(1)*x(2);   

dx(5) = alpha*(x(2)^2); 


%(semi VOID)