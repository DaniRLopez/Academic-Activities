%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 1.4.3 Numerical Solution of Rayleigh–Lamb Frequency Equations

Abstract
This file contains the solution of the eq. 1.92 for the term LHSA using the
algorithm presented in pg. 27 that implements the bisection method for the
solution of the eq. 1.91 for antisymmetric modes.

D.R Lopez
2020

%}

%% Cleaning the workspace
clc
clear
close all
hold on
%% Plotting the Phase velocity figure for antisymmetric modes
error = 1e-10;                        % tolerance for the bisection method
cs = 3.2e3;                           % shear wave velocity 
cl = 6.3e3;                           % longitudinal wave velocity
cp  = 1:100:10e3;                     % c points
WD  = 2*pi()*[1e+0:1e+1:20e+3];       % frequency range

for u = 1:length(WD)
wd = WD(u);

% Eq. 1.92 that represents the modified symbols p and q presented in eq.
% 1.66
LHSA = zeros(length(cp),1);
roots = [];
g     = 1;

for k = 1:length(cp)
    c = cp(k);
    qg = sqrt((1/cs^2)-(1/c^2));       %     eq. 1.93
    pg = sqrt((1/cl^2)-(1/c^2));       %     eq. 1.93
    % Eq. 1.92 for antisymmetric modes
    term1 = qg*tan(qg*wd);
    term2 = ((qg^2-(1/c^2))^2)*tan(pg*wd);
    term3 = 4*pg*c^2;
    LHSA(k)  = term1 + (term2/term3);  %     eq. 1.92
    
    if k ~= 1
       %% Bisection method
      if sign(LHSA(k)) ~= sign(LHSA(k-1))
           a = cp(k-1);
           b = cp(k);
           tol = b-a;

           while tol > error
           % d is the new point that represents the possible root of the
           % function, in this case,LHSA
           d = (a+b)/2;
           % function evaluated in the points a,b and d
           qg1 = sqrt((1/cs^2)-(1/a^2));   
           pg1 = sqrt((1/cl^2)-(1/a^2));
           term11 = qg1*tan(qg1*wd);
           term21 = ((qg1^2-(1/a^2))^2)*tan(pg1*wd);
           term31 = 4*pg1*a^2;
           f1  = term11 + (term21/term31);
           
           qg2 = sqrt((1/cs^2)-(1/b^2));   
           pg2 = sqrt((1/cl^2)-(1/b^2));
           term12 = qg2*tan(qg2*wd);
           term22 = ((qg2^2-(1/b^2))^2)*tan(pg2*wd);
           term32 = 4*pg2*b^2;
           f2  = term12 + (term22/term32);
           
           qg3 = sqrt((1/cs^2)-(1/d^2));   
           pg3 = sqrt((1/cl^2)-(1/d^2));
           term13 = qg3*tan(qg3*wd);
           term23 = ((qg3^2-(1/d^2))^2)*tan(pg3*wd);
           term33 = 4*pg3*d^2;
           f3  = term13 + (term23/term33);
           % Criteria for bisection method
               if sign(f1) == sign(f3)
                   a = d;
               else
                   b = d;
               end
           % difference between new a and b
           tol = b-a;
           end
           roots(g,:) = [wd,d];
           g = g+1;
      else
          
      end
    end
end
figure(1)
hold on
if length(roots()) > 0
%% Figure 1.12 for antisymmetric waves
plot(roots(:,1)/(2*pi()*(1e3)),roots(:,2)/(1e3),'.b')
end
end

grid on
set(gca,'TickLabelInterpreter','latex')
xlabel('Fequency parameter $f \cdot d [MHz \cdot m]$','Interpreter','latex')
ylabel('Phase velocity $C_p$','Interpreter','latex')

%% Bye! 
return