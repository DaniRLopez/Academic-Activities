%{
Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019
Section : 2.4.2 Gauss Quadrature
Abstract
This file contains values of weights and abscissae of Gauss quadrature
using the Legendre Polynomials. This program replicates table 2.9
until n = 6. If you want to obtain another n with their respective 
abscissae and weights only must change the line number 26 with your desired 
value n plus number 1.
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
ninitial = 0;
nfinal   = 7;
n  = ninitial:1:nfinal;
Pn = sym(zeros(1,length(n)));
for x = 1:length(n)
    first  = 1/((2^n(x))*(factorial(n(x))));
    second = diff(((xi^2)-1)^n(x),n(x));
    Pn(x)  = first*second;
end
%% Abscissae for Gauss Quadratures
%  The abscissae used for to calculate the weigths of Gauss Quadrature
%  are the roots of the Legendre Polynomials with order q, according with
%  eq. 2.55
Abscissae = sym(zeros(length(n)-2,length(n)-2));
for w = 1:(length(n)-2)
    eqn    = Pn(w+2);
    Abscissae(w,1:w+1) = vpa(solve(eqn,xi)).';
end
%% Weights for Gauss Quadratures
% The weights are calculated using  the eq. 2.54 for n = 1 to 6.
wG = sym(zeros(length(n)-1,length(n)-1));
fprintf('-Values of Abscissae and Weights for Gauss Quadratures-\n');
fprintf('Element order   Abscissae coordinates            Weights\n');
for k = 1:(length(n)-2)
    s = k+2;
    fprintf('n = %d',k);
    for g = 1:(s-1)
        second = (1-(Abscissae(k,g))^2)*((subs(diff(Pn(s)),Abscissae(k,g)))^2);
        wG(k,g) = 2/second;
        if g == 1 && Abscissae(k,g) >= 0
        formatSpect = '            a%d = %4.6f                w%d = %4.6f\n';
        fprintf(formatSpect,[g,Abscissae(k,1),g,wG(k,1)]);
        elseif g == 1 && Abscissae(k,g) < 0
        formatSpect = '            a%d = %4.6f               w%d = %4.6f\n';
        fprintf(formatSpect,[g,Abscissae(k,1),g,wG(k,1)]);
        elseif Abscissae(k,g) < 0
        formatSpecta = '                 a%d = %4.6f               w%d = %4.6f\n';
        fprintf(formatSpecta,[g,Abscissae(k,g),g,wG(k,g)]);   
        else
        formatSpect1 = '                 a%d = %4.6f                w%d = %4.6f\n';
        fprintf(formatSpect1,[g,Abscissae(k,g),g,wG(k,g)]);
        end
    end
end
%% End
return