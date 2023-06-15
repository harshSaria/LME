set(groot,'defaultLineLineWidth',2.0)

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
n=4;
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

vary_d1=d1_vary(linspace(2,9,20)*10^(-3),d2,l1,l2,mdot_lox,rho_lox,mu_lox, a_rough);
% vary_d2=d2_vary(d1,linspace(2,6,20)*10^(-3),l1,l2,mdot_lox,rho_lox,mu_lox, a_rough);

function f = friction_factor(Re, r_rough)
f=(-2.*log10((r_rough./3.7)-(5.02./Re).*log10(r_rough-(5.02./Re).*((r_rough./3.7)+(13./Re))))).^(-2);
% f=(-1.8.*log10((r_rough./3.7).^1.11 + (6.9./Re))).^(-2);
% f=(1.14-2.*log10(r_rough+(21.25./Re.^0.9))).^(-2);
% f=(-2.*log10((r_rough./3.7)+(5.74./Re.^0.9))).^(-2);
end

function vary1 = d1_vary(d1,d2,l1,l2,mdot_lox,rho_lox,mu_lox, a_rough)
for i=1:length(mdot_lox)
a1=3.14.*d1.*d1;
v1=mdot_lox(i)./(a1*rho_lox);
r_rough=0.5*a_rough./d1;
Re1=rho_lox*v1*2.*d1/mu_lox;
f1=friction_factor(Re1,r_rough);

delta_P1=f1.*l1.*v1.*v1./(2*9.81*2.*d1);

figure(1)
grid on
xlabel("Radius (m)")
ylabel("Pressure Drop (bar)")
title("Pressure Drop vs Radius")
plot(d1,delta_P1)
scatter(d1,delta_P1)
hold on

figure(2)
grid on
xlabel("Radius (m)")
ylabel("Friction Factor")
title("Friction Factor vs Radius")
plot(d1,f1)
scatter(d1,f1)
hold on

figure(3)
grid on
xlabel("Reynolds Number")
ylabel("Friction Factor")
title("Reynolds Number vs Friction Factor")
plot(Re1,f1)
scatter(Re1,f1)
hold on

figure(4)
grid on
xlabel("Radius (m)")
ylabel("Relative Roughness (m^{-1})")
title("Relative Roughness vs Radius")
plot(d1,r_rough)
scatter(d1,r_rough)
hold on

a2=3.14*d2*d2;
v2=mdot_lox(i)/(a2*rho_lox);
r_rough=0.5*a_rough/d2;
Re2=rho_lox*v2*2*d2/mu_lox;
f2=friction_factor(Re2,r_rough);

delta_P2=f2.*l2.*v2.*v2./(2*9.81*2*d2);
end
end

% function vary2 = d2_vary(d1,d2,l1,l2,mdot_lox,rho_lox,mu_lox, a_rough)
% for i=1:length(mdot_lox)
% a2=3.14.*d2.*d2;
% v2=mdot_lox(i)./(a2*rho_lox);
% r_rough=0.5*a_rough./d2;
% Re2=rho_lox*v2*2.*d2/mu_lox;
% f2=friction_factor(Re2,r_rough);
% 
% delta_P2=f2.*l2.*v2.*v2./(2*9.81*2.*d2)
% 
% figure(1)
% grid on
% xlabel("Radius (m)")
% ylabel("Pressure Drop (bar)")
% title("Pressure Drop vs Radius")
% scatter(d2,delta_P2)
% plot(d2,delta_P2)
% hold on
% 
% figure(2)
% grid on
% xlabel("Radius (m)")
% ylabel("Friction Factor")
% title("Friction Factor vs Radius")
% scatter(d2,f2)
% plot(d2,f2)
% hold on
% 
% figure(3)
% grid on
% xlabel("Reynolds Number")
% ylabel("Friction Factor")
% title("Reynolds Number vs Friction Factor")
% scatter(Re2,f2)
% plot(Re2,f2)
% hold on
% 
% figure(4)
% grid on
% xlabel("Radius (m)")
% ylabel("Relative Roughness (m^{-1})")
% title("Relative Roughness vs Radius")
% scatter(d2,r_rough)
% plot(d2,r_rough)
% hold on
% 
% a1=3.14*d1*d1;
% v1=mdot_lox(i)/(a1*rho_lox);
% r_rough=0.5*a_rough/d1
% Re1=rho_lox*v1*2*d1/mu_lox;
% f1=friction_factor(Re1,r_rough);
% 
% delta_P1=f1.*l1.*v1.*v1./(2*9.81*2*d1);
% end
% end
