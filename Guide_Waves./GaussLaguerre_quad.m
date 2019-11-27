%{
Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019
Section : 2.4.3 Gauss-Laguerre Quadrature
Abstract
This file contains values of weights and Gnc of Gauss quadrature
using the full Chebyshev Polynomials. This program replicates table 2.9
until n = 6. If you want to obtain another n with their respective 
Gnc and weights only must change the line number 26 with your desired 
value n plus number 2.
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
n  = 0:1:8;
Gn = sym(zeros(1,length(n)));
for x = 1:length(n)
    first  = exp(xi)/factorial(n(x));
    second = diff((exp(-xi))*(xi^n(x)),n(x));
    Gn(x)  = first*second;
end
%% Gnc for Gauss-Laguerre Quadratures
%  The Gnc used for to calculate the weigths of Gauss Quadrature
%  are the roots of the Laguerre Polynomials with order q, according
%  with eq. 2.59
Gnc = sym(zeros(length(n)-2,length(n)-2));
for w = 1:(length(n)-2)
    eqn    = (Gn(w+2));
    Gnc(w,1:w+1) = vpa(solve(eqn,xi)).';
end
%% Weights for Gauss-Laguerre Quadratures
% The weights are calculated using  the eq. 2.54 for n = 1 to 6.
wGL = sym(zeros(length(n)-1,length(n)-2));
fprintf('-Values of Gnc and Weights for Gauss-Laguerre Quadratures-\n');
fprintf('Element order   Gnc coordinates                 Weights\n');
for k = 1:(length(n)-3)
    s = k+2;
    fprintf('n = %d',k);
    for g = 1:(s-1)
        second = ((s)^2)*((subs(Gn(s+1),Gnc(k,g)))^2);
        wGL(k,g) = Gnc(k,g)/second;
        if g == 1
        formatSpect = '            a%d = %4.6f               w%d = %4.6f\n';
        fprintf(formatSpect,[g,Gnc(k,1),g,wGL(k,1)]);
        elseif Gnc(k,g) >10 
        formatSpecta = '                 a%d = %4.6f              w%d = %4.6f\n';
        fprintf(formatSpecta,[g,Gnc(k,g),g,wGL(k,g)]);   
        else
        formatSpect1 = '                 a%d = %4.6f               w%d = %4.6f\n';
        fprintf(formatSpect1,[g,Gnc(k,g),g,wGL(k,g)]);
        end
    end
end
%% End
return