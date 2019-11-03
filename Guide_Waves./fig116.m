%{

Book    : Guided Waves in Structures for SHM - the time-domain 
          spectral element method
Authors : Wieslaw Ostachowicz; Pawel Kudela; Marek Krawczuk; Arkadiusz Zak
Year    : 2012
Matlab  : 2019

Section : 1.4.5 Shear Horizontal Waves.


Abstract
This file contains the solutions of phase velocity dispersion curves for 
shear horizontal (SH) waves, presented in eq.1.108 for both symmetric and 
antisymmetric modes, presented in fig 1.16.


D.R. Lopez
2019


%}

%% File settings
clc
clear
close all

%% Graphics settings
hold on
ylim([0 10]);
xlabel('Fequency parameter $f \cdot d [MHz \cdot m]$','Interpreter','latex')
ylabel('Phase velocity $C_p$','Interpreter','latex')

%% Explicit solutions

n = 0:1:5;

eta_S_d = n*pi();                     % solution eq. 105
eta_A_d = (2*n+1)*pi()/2;             % solution eq. 106
cs      = 3200;                       % shear wave velocity
fd      = (10^3)*(0.01:0.001:20);     % Frecuency parameter
wd      = pi()*fd;
T       = zeros(length(wd),length(eta_S_d));

%% Phase velocity for symmetric SH waves
%  In this part is solved the eq. 1.108 for
%  symmetric modes. The symmetric modes of SH waves
%  are obtained by solving the eq. 1.105. For a shorter
%  computing time some restricctions are imposed.

Cp = zeros(length(wd),length(eta_S_d));
for y = 1:length(wd)
    for g = 1:length(eta_S_d)
        first = sqrt(1-((eta_S_d(g)^2)*(cs/wd(y))^2));% eq 1.108
        Cp(y,g) = cs/(first*1000);
        T(y,g)  = fd(y)/(10^3);
        if wd(y) == cs || isreal(Cp(y,g)) == 0
            Cp(y,g) = 0;
            T(y,g)  = 0;
        end
    end
end

%% Phase velocity for antisymmetric SH waves
%  In this part is solved the eq. 1.108 for
%  antisymmetric modes. The antisymmetric modes of SH waves
%  are obtained by solving the eq. 1.106. For a shorter
%  computing time some restricctions are imposed.

Cpa = zeros(length(wd),length(eta_A_d));
Ta  = zeros(length(wd),length(eta_A_d));
for y = 1:length(wd)
    for g = 1:length(eta_A_d)
        first = sqrt(1-((eta_A_d(g))^2*(cs/wd(y))^2));% eq 1.108
        Cpa(y,g) = cs/(first*1000);
        Ta(y,g)  = fd(y)/(10^3);
        if wd(y) == cs || isreal(Cpa(y,g)) == 0
            Cpa(y,g) = 0;
            Ta(y,g)  = 0;
        end
    end
end

%% Phase velocity graphics 
q = plot(T(:,1),Cp(:,1),'.b','MarkerSize',2);
plot(T(:,2),Cp(:,2),'.b',T(:,3),Cp(:,3),'.b',T(:,4),...
    Cp(:,4),'.b',T(:,5),Cp(:,5),'.b','MarkerSize',2);

r = plot(Ta(:,1),Cpa(:,1),'.r','MarkerSize',2);
plot(Ta(:,2),Cpa(:,2),'.r',Ta(:,3),Cpa(:,3),'.r',Ta(:,4),...
    Cpa(:,4),'.r',Ta(:,5),Cpa(:,5),'.r','MarkerSize',2);

legend([q r],{'Symmetric','Antisymmetric'})

grid on

%%
return;
