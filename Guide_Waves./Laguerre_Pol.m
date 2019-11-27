%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 2.1.3 Laguerre Polynomials


Abstract
This file contains the original form of Laguerre Polynomial of the nth order
usign de Rodrigues's formula for to generate the polynomials. Also, in this
script are calculated the full Laguerre Polynomial with their roots in
the normalised coordinated system [0, +Inf].

In this program are calculated the Laguerre Polynomials of 0th to 5th grade
and n = 1 to n = 6 Full Laguerre Polynomials.

D.R. Lopez
2019


%}
%% 
clc
clear
close all
%% Loading the syms toolbox
syms xi 
%% Laguerre Formula
%  In this section is used the eq. 2.15 for to calculate the Laguerre
%  Polynomial
n  = 0:1:7;
Gn = sym(zeros(1,length(n)));
for x = 1:length(n)
    first  = exp(xi)/factorial(n(x));
    second = diff((exp(-xi))*(xi^n(x)),n(x));
    Gn(x)  = first*second;
end
%% Graphs of the Laguerre Polynomials
%  This section replicates fig. 2.6
figure(1)
hold on
axis([0 20 -10 20])
for z = 1:(length(n))-2
    fplot(Gn(z),[0,20])
end
%  set plot in LaTeX format
set(gca,'TickLabelInterpreter','latex');
legend('n=0','n=1','n=2','n=3','n=4','n=5','Location','Best',...
    'Interpreter','latex');
xlabel('$\xi$- coordinate','Interpreter','latex')
ylabel('Laguerre Polynomial $G_n(\xi)$','Interpreter','latex')
grid on
%% Complete Laguerre Polynomials
%  This section replicates the table 2.6 with the roots of the eq. 2.18
Gnc = sym(zeros(length(n)-2,length(n)-2));
fprintf('-Values of Laguerre nodes in the normalised coordinate system-\n');
fprintf('Element order     Node coordinates\n');
for w = 1:(length(n)-2)
    eqn    = (xi)*(Gn(w+1));
    Gnc(w,1:w+1) = vpa(solve(eqn,xi)).';
    fprintf('n = %d',w);
    fprintf('              xi = %3f \n',Gnc(w,1));
    fprintf('                   xi = %3f \n',Gnc(w,2:w+1));
    fprintf('\n');
end
return;
