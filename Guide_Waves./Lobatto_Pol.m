%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 2.1.1 Lobatto Polynomials


Abstract
This file contains the original form of Lobatto Polynomial of the nth order
usign de Rodrigues's formula for to generate the polynomials. Also, in this
script are calculated the complete Lobatto Polynomial with their roots in
the normalised coordinated system [-1,+1].

In this program are calculated the Lobatto Polinomials of 0th to 5th grade.

D.R. Lopez
2019


%}

%% 
clc
clear
close all

%% Loading the syms toolbox
syms xi 

%% Rodrigues's formula
%  In this section is used the eq. 2.5 for to calculate the Legendre
%  Polynomial
n  = 0:1:6;
Pn = sym(zeros(1,length(n)));
for x = 1:length(n)
    first  = 1/((2^n(x))*(factorial(n(x))));
    second = diff(((xi^2)-1)^n(x),n(x));
    Pn(x)  = first*second;
end
%% Lobatto Polynomials
%  In this section is used the eq. 2.6 for to calculate the Legendre
%  Polynomial. This algorithm replicates Table 2.1
Ln = sym(zeros(1,length(n)-1));
for y = 1:(length(n)-1)
    Ln(y) = diff(Pn(y+1),xi);
end
%% Graphs of the Lobatto Polynomials
%  This section replicates fig. 2.4
figure(1)
hold on
for z = 1:(length(n)-1)
    fplot(Ln(z),[-1,+1])
end
%  set plot in LaTeX format
set(gca,'TickLabelInterpreter','latex');
legend('n=0','n=1','n=2','n=3','n=4','n=5','Location','Best',...
    'Interpreter','latex');
xlabel('$\xi$- coordinate','Interpreter','latex')
ylabel('Lobatto Polynomial $L_n(\xi)$','Interpreter','latex')
grid on
%% Complete Lobatto Polynomials
%  This section replicates the table 2.2 with the roots of the eq. 2.9
Lnc = sym(zeros(length(n)-1,length(n)-1));
fprintf('-Values of Lobatto nodes in the normalised coordinate system-\n');
fprintf('Element order     Node coordinates\n');
for w = 1:(length(n)-1)
    eqn    = (1-(xi^2))*(Ln(w));
    Lnc(w,1:w+1) = vpa(solve(eqn,xi)).';
    fprintf('n = %d',w);
    fprintf('              xi = %3f \n',Lnc(w,1));
    fprintf('                   xi = %3f \n',Lnc(w,2:w+1));
    fprintf('\n');
end
return;