%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 2.1.2 Chebyshev Polynomials


Abstract
This file contains the original form of Chebyshev Polynomial of the of the
second type using the recursive formula for to generate the polynomials.
Also, in this script are calculated the full Chebyshev Polynomial with 
their roots in the normalised coordinated system [-1,+1].

In this program are calculated the Chebyshev Polinomials of 0th to 5th
grade.

D.R. Lopez
2019
%}

%% 
clc
clear
close all
%% Loading the syms toolbox
syms xi 
%% Chebyshev Polynomials
%  In this section is used the eq. 2.10 for to calculate the recursive
%  formunal. This algorithm replicates Table 2.3
n = 0:1:5;
Un = sym(zeros(1,length(n)));
for y = 1:(length(n))
    if n(y) == 0
        Un(1) = 1;
    elseif n(y) == 1
        Un(2) = 2*xi;
    else
        Un(y) = (2*xi*Un(y-1))-Un(y-2);
    end
end
%% Graphs of the Chebyshev Polynomials
%  This section replicates fig. 2.5
figure(1)
hold on
for z = 1:(length(n))
    fplot(Un(z),[-1,+1])
end
%  set plot in LaTeX format
set(gca,'TickLabelInterpreter','latex');
legend('n=0','n=1','n=2','n=3','n=4','n=5','Location','Best',...
    'Interpreter','latex');
xlabel('$\xi$- coordinate','Interpreter','latex')
ylabel('Chebyshev Polynomial $U_n(\xi)$','Interpreter','latex')
grid on
%% Complete Lobatto Polynomials
%  This section replicates the table 2.4 with the roots of the eq. 2.13
Tnc = sym(zeros(length(n),length(n)));
fprintf('-Values of Chebyshev nodes in the normalised coordinate system-\n');
fprintf('Element order     Node coordinates\n');
for w = 1:(length(n))
    eqn    = (1-(xi^2))*(Un(w));
    Tnc(w,1:w+1) = vpa(solve(eqn,xi)).';
    fprintf('n = %d',w);
    fprintf('              xi = %3f \n',Tnc(w,1));
    fprintf('                   xi = %3f \n',Tnc(w,2:w+1));
    fprintf('\n');
end
%%
return