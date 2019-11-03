%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 1.4.3 Numerical solutions of Rayleigh-Lamb Frequency Equations


Abstract
This file contains the solutions of phase velocity dispersion curves for 
symmetric modes presented in eq 1.91. For a shorter computation time
some resctrictions are imposed in the proposed method.

Convergence criterion for the Bisection Method: 1e-7

D.R. Lopez
2019


%}

%% 
clc
clear
close all

%% Graphics settings
hold on
xlabel('Fequency parameter $f \cdot d  [MHz \cdot m]$','Interpreter','latex')
ylabel('Phase velocity $C_p$','Interpreter','latex')
ylim([0 10])

cs = 3.2e3;                                  % shear wave velocity 
cl = 6.3e3;                                  % longitudinal wave velocity
c_sup  = 0.01:1:10e3;

frecd = 10e3*(0.0001:0.05:2.0);
C = zeros(1,2);
D = zeros(size(frecd,2),size(c_sup,2));

%% LHS equation solution

% in this section is solved the eq. 1.92 for antisymmetric modes using the
% bisection method proposed in pg. 27. For a better computing time some 
% restrictions are imposed.

for x = 1:size(frecd,2)
fd = frecd(x);
wd = fd*pi();
h = 1;
u = 1;
for y = 1:size(c_sup,2)
   h  = 1; 
   u  = 1;
   c = c_sup(y);
 
while h == 1

% Eq. 1.92 that represent the modified symbols p and q presented in eq.
% 1.66
qg = sqrt((1/cs^2)-(1/c^2));   
pg = sqrt((1/cl^2)-(1/c^2));

% Eq. 1.92 for symmetric modes
term1 = 4*pg*tan(pg*wd);
term2 = ((qg^2-(1/c^2))^2)*(c^2);
term3 = tan(qg*wd);
LHS  = (term3/qg) + (term1/term2);

   if c > 10e3
       % D(x,y) = c;
       break
   end

if u == 1
    C(1) = c;
    u = 2;
    eas(1) = LHS;
    c = c + c/1000;
else
    if sign(LHS) == sign(eas(1))
        eas(1) = LHS;
        C(1) = c;
        c = c + c/1000;
        C(2) = c;
    else
        % Bisection method in pg. 27
        while abs(C(1)-C(2)) > 1e-7
            if eas(1) < LHS
                C(1) = C(1);
                C(2) = (C(2)+C(1))/2;
                c = C(2);
            elseif eas(1) > LHS
                C(2) = C(2);
                C(1) = (C(2)+C(1))/2;
                c = C(1);
            end
            
            qg = sqrt((1/cs^2)-(1/c^2));
            pg = sqrt((1/cl^2)-(1/c^2));

            term1 = 4*pg*tan(pg*wd);
            term2 = ((qg^2-(1/c^2))^2)*(c^2);
            term3 = tan(qg*wd);
            LHS  = (term3/qg) + (term1/term2);
        
        end
        
        if abs(C(1)-C(2))<1e-7
            D(x,y) = C(2);
            h = 2;
        end

    end
end
end

end
hold on
plot(fd*(ones(1,size(D,2)))/1000,D(x,:)/1000,'.b','MarkerSize',2)
end

%%
return;

