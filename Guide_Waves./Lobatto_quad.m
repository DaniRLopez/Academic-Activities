%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 2.4.1 Lobatto Quadrature


Abstract
This file contains values of weights and abscissae of Lobatto quadrature
using the full Lobatto Polynomials. This program replicates table 2.8 until
n = 6, if you want to obtain another n with their respective abscissae and
weights only must change the line number 34 with your desired value.


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
nfinal   = 6;
n  = ninitial:1:nfinal;
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
%% Complete Lobatto Polynomials
%  This section replicates the table 2.2 with the roots of the eq. 2.9
Lnc = sym(zeros(length(n)-1,length(n)-1));
for w = 1:(length(n)-1)
    eqn    = (1-(xi^2))*(Ln(w));
    Lnc(w,1:w+1) = vpa(solve(eqn,xi)).';
end
%% Abscissae for Lobatto Quadratures
%  The abscissae used for to calculate the weigths of Lobatto Quadrature
%  are the same points found for the complete lobatto Polynomials. In this
%  way the abscissae used in this program will be the xi for complete
%  (full) Lobatto Polynomials.
abscissae = Lnc;
%% Weights for Lobatto Quadratures
% The weights are calculated using  the eq. 2.50 for n = 1 to 6.
wL = sym(zeros(length(n)-1,length(n)-1));
fprintf('-Values of Abscissae and Weights for Lobatto Quadratures-\n');
fprintf('Element order   Abscissae coordinates      Weights\n');
for k = 1:(length(n)-1)
    s = k+1;
    fprintf('n = %d',k);
    for g = 1:s
        second = s*(s-1)*((subs(Pn(s),Lnc(k,g)))^2);
        wL(k,g) = 2/second;
        if g == 1
        formatSpect = '            a%d = %4.6f               w%d = %4.6f\n';
        fprintf(formatSpect,[g,Lnc(k,1),g,wL(k,1)]);
        elseif Lnc(k,g) < 0
        formatSpecta = '                 a%d = %4.6f               w%d = %4.6f\n';
        fprintf(formatSpecta,[g,Lnc(k,g),g,wL(k,g)]);   
        else
        formatSpect1 = '                 a%d = %4.6f                w%d = %4.6f\n';
        fprintf(formatSpect1,[g,Lnc(k,g),g,wL(k,g)]);
        end
    end
end
%% End
return