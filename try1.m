% knowns
thrust_level=[55,100,115];
mdot_lox=[0.4044,0.7317,0.841];
mdot_ch4=[0.1072,0.1939,0.223];
mr=[3.773,3.773,3.773];
t_lox=92;
t_ch4=252;
p_lox=158;
p_ch4=145;
rho_lox=1164.8;
mu_lox=215.39e-6;
rho_ch4=2;
mu_ch4=3;

% material property
a_rough = 0.03e-3;

% geometric parameters (m)
hole_dia=2e-3;
n=4e-3;
d1=4.5e-3;
d2=3e-3;
l1=14e-3;
l2=36.2e-3;
lip=0.9e-3;
e=3.4e-3;
df1=0.4e-3;
df2=0.8e-3;
lf1=9.7e-3;
lf2=9.8e-3;


a1=3.14*d1*d1;
v1=mdot_lox/(a1*rho_lox);
r_rough=0.5*a_rough/d1;
Re1=rho_lox*v1*2*d1/mu_lox;
f1=friction_factor(Re1,r_rough);

delta_P1=f1.*l1.*v1.*v1./(2*9.81*2*d1);
cd1=mdot_lox./(a1.*sqrt(2*rho_lox*delta_P1.*0.1.*10^6));

a2=3.14*d2*d2;
v2=mdot_lox/(a2*rho_lox);
r_rough=0.5*a_rough/d2;
Re2=rho_lox*v2*2*d2/mu_lox;
f2=friction_factor(Re2,r_rough);

delta_P2=f2.*l2.*v2.*v2./(2*9.81*2*d2);
cd2=mdot_lox./(a2*sqrt(2*rho_lox.*delta_P2.*0.1.*10^6));

total_delP_lox=delta_P2+delta_P1; % using friction factor formula
total_cd_lox=cd1+cd2;

% using simple formula mdot=cd*A*sqrt(2*rho*delP)
cdguess=1;
delp1=((mdot_lox./a1).^2)./(cdguess*2*rho_lox*10^5);
delp2=((mdot_lox./a2).^2)./(cdguess*2*rho_lox*10^5);
tot_delP = delp1+delp2;

diff = total_cd_lox-tot_delP; % difference between delta P value from two method

function f = friction_factor(Re, r_rough)
f=(-2.*log10((r_rough./3.7)-(5.02./Re).*log10(r_rough-(5.02./Re).*((r_rough./3.7)+(13./Re))))).^(-2);
% f=(-1.8.*log10((r_rough./3.7).^1.11 + (6.9./Re))).^(-2);
% f=(1.14-2.*log10(r_rough+(21.25./Re.^0.9))).^(-2);
% f=(-2.*log10((r_rough./3.7)+(5.74./Re.^0.9))).^(-2);
end


