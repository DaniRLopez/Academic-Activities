%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 1.4.3 Numerical Solution of Rayleigh–Lamb Frequency Equations

Abstract
This file contains the solution of the eq. 1.92 for the term LHSS using the
algorithm presented in pg. 27 that implements the bisection method for the
solution of the eq. 1.91 for symmetric modes.

D.R Lopez
2020

%}

%% Cleaning the workspace
clc
clear
close all
hold on
%% Plotting the Phase velocity figure for symmetric modes
error = 1e-10;                        % tolerance for the bisection method
cs = 3.2e3;                           % shear wave velocity 
cl = 6.3e3;                           % longitudinal wave velocity
cp  = 1:100:10e3;                     % c points
WD  = 2*pi()*[1e+0:1e+1:20e+3];       % frequency range

for u = 1:length(WD)
wd = WD(u);

% Eq. 1.92 that represents the modified symbols p and q presented in eq.
% 1.66

LHSS = zeros(length(cp),1);
roots = [];
g     = 1;
% Eq. 1.92 for symmetric modes
for k = 1:length(cp)
    c = cp(k);
    qg = sqrt((1/cs^2)-(1/c^2));            %    eq. 1.93
    pg = sqrt((1/cl^2)-(1/c^2));            %    eq. 1.93
    % Eq. 1.92 for symmetric modes
    term1 = tan(qg*wd)/qg;
    term2 = 4*pg*tan(pg*wd);
    term3 = ((qg^2-(1/c^2))^2)*(c^2);
    LHSS(k)  = term1 + (term2/term3);       %    eq. 1.92
    
    if k ~= 1
      if sign(LHSS(k)) ~= sign(LHSS(k-1))
          %% Bisection method
           a = cp(k-1);
           b = cp(k);
           tol = b-a;
           while tol > error
           % d is the new point that represents the possible root of the
           % function, in this case,LHSS
           d = (a+b)/2;
           % function evaluated in the points a,b and d
           qg1 = sqrt((1/cs^2)-(1/a^2));   
           pg1 = sqrt((1/cl^2)-(1/a^2));
           term11 = tan(qg1*wd)/qg1;
            term21 = 4*pg1*tan(pg1*wd);
           term31 = ((qg1^2-(1/a^2))^2)*(a^2);
           f1  = term11 + (term21/term31);
           
           qg2 = sqrt((1/cs^2)-(1/b^2));   
           pg2 = sqrt((1/cl^2)-(1/b^2));
           term12 = tan(qg2*wd)/qg2;
            term22 = 4*pg2*tan(pg2*wd);
           term32 = ((qg2^2-(1/b^2))^2)*(b^2);
           f2  = term12 + (term22/term32);
           
           qg3 = sqrt((1/cs^2)-(1/d^2));   
           pg3 = sqrt((1/cl^2)-(1/d^2));
           term13 = tan(qg3*wd)/qg3;
            term23 = 4*pg3*tan(pg3*wd);
           term33 = ((qg3^2-(1/d^2))^2)*(d^2);
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
plot(roots(:,1)/(2*pi()),roots(:,2),'.b')
end
end

grid on
set(gca,'TickLabelInterpreter','latex')
xlabel('Fequency parameter $f \cdot d [MHz \cdot m]$','Interpreter','latex')
ylabel('Phase velocity $C_p$','Interpreter','latex')

%% Bye!
return
