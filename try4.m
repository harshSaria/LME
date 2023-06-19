% d1, d2, l1 and l2 are paramterized for axial lox flow

set(groot,'defaultLineLineWidth',2.0)

% user input
vary = input("What do you want to vary? Type 1 for d1, 2 for d2, 3 for L1, 4 for L2: ");
start=input("Enter the starting value of parameter you want to vary in metres: ");
final=input("Enter the ending value of parameter you want to vary in metres: ");
data_points=input("Enter number of data points you want between them: ");

% knowns with default values
thrust_level=[55,100,115];
mdot_lox=[0.4044,0.7317,0.841];
mdot_ch4=[0.1072,0.1939,0.223];
mr=[3.773,3.773,3.773];
t_lox=linspace(92,92,data_points);
t_ch4=linspace(252,252,data_points);
p_lox=linspace(158,158,data_points);
p_ch4=linspace(145,145,data_points);
rho_lox=linspace(1164.8,1164.8,data_points);
mu_lox=linspace(215.39e-6, 215.39e-6, data_points);
rho_ch4=linspace(2,2,data_points);
mu_ch4=linspace(3,3,data_points);

% material property
a_rough = linspace(0.03e-3,0.03e-3,data_points);

% geometric parameters with default values (m)
hole_dia=linspace(2e-3,2e-3,data_points);
n=linspace(4,4,data_points);
d1=linspace(4.5e-3,4.5e-3,data_points);
d2=linspace(3e-3,3e-3,data_points);
l1=linspace(14e-3,14e-3,data_points);
l2=linspace(36.2e-3,36.2e-3,data_points);
lip=linspace(0.9e-3,0.9e-3,data_points);
e=linspace(3.4e-3,3.4e-3,data_points);
df1=linspace(0.4e-3,0.4e-3,data_points);
df2=linspace(0.8e-3,0.8e-3,data_points);
lf1=linspace(9.7e-3,9.7e-3,data_points);
lf2=linspace(9.8e-3,9.8e-3,data_points);

% choice of variation
if vary==1
    d1=linspace(start,final,data_points);
    variable=d1;
elseif vary==2
    d2=linspace(start, final, data_points);
    variable=d2;
elseif vary==3
    l1=linspace(start, final, data_points);
    variable=l1;
elseif vary==4
    l2=linspace(start, final, data_points);
    variable=l2; 
end

% calling
vary_d1=d1_vary(d1,d2,l1,l2,mdot_lox,rho_lox,mu_lox, a_rough, vary,variable);

function f = friction_factor(Re, r_rough)
f=(-2.*log10((r_rough./3.7)-(5.02./Re).*log10(r_rough-(5.02./Re).*((r_rough./3.7)+(13./Re))))).^(-2);
% f=(-1.8.*log10((r_rough./3.7).^1.11 + (6.9./Re))).^(-2);
% f=(1.14-2.*log10(r_rough+(21.25./Re.^0.9))).^(-2);
% f=(-2.*log10((r_rough./3.7)+(5.74./Re.^0.9))).^(-2);
end

function vary1 = d1_vary(d1,d2,l1,l2,mdot_lox,rho_lox,mu_lox, a_rough,vary, variable)
for i=1:length(mdot_lox)
% first pipe
a1=3.14.*d1.*d1;
v1=mdot_lox(i)./(a1.*rho_lox);
r_rough1=0.5.*a_rough./d1;
Re1=rho_lox.*v1.*2.*d1./mu_lox;
f1=friction_factor(Re1,r_rough1);
delta_P1=f1.*l1.*v1.*v1./(2*9.81*2.*d1);

% second pipe
a2=3.14.*d2.*d2;
v2=mdot_lox(i)./(a2.*rho_lox);
r_rough2=0.5.*a_rough./d2;
Re2=rho_lox.*v2.*2.*d2./mu_lox;
f2=friction_factor(Re2,r_rough2);
delta_P2=f2.*l2.*v2.*v2./(2*9.81*2.*d2);

% combining both pipe data
delta_P12=[delta_P1',delta_P2'];
f12=[f1',f2'];
Re12=[Re1',Re2'];
r_rough12=[r_rough1',r_rough2'];

delta_P_total = delta_P2 + delta_P1; % total pressure drop summing of two pipes

% Non Dimensional Number
dia_rat = d1./d2;
lengh_dia = l2./d2;

% plotting
if vary==1 || vary==3 % that is for pipe 1
figure(1)
grid on
xlabel("Variable (m)")
ylabel("Pressure Drop (bar)")
title("Pressure Drop vs Variable")
scatter(variable,delta_P12(:,1))
plot(variable,delta_P12(:,1))
hold on


figure(2)
grid on
xlabel("Variable (m)")
ylabel("Friction Factor")
title("Friction Factor vs Variable")
plot(variable,f12(:,1))
scatter(variable,f12(:,1))
hold on


figure(3)
grid on
xlabel("Reynolds Number")
ylabel("Friction Factor")
title("Reynolds Number vs Friction Factor")
plot(Re12(:,1),f12(:,1))
scatter(Re12(:,1),f12(:,1))
hold on


figure(4)
grid on
xlabel("Variable (m)")
ylabel("Relative Roughness (m^{-1})")
title("Relative Roughness vs Variable")
plot(variable,r_rough12(:,1))
scatter(variable,r_rough12(:,1))
hold on


figure(5)
grid on
xlabel("Variable (m)")
ylabel("Reynolds Number")
title("Reynolds Number vs Variable")
plot(variable,Re12(:,1))
scatter(variable,Re12(:,1))
hold on

figure(6)
grid on
xlabel("D_1/D_2")
ylabel("Pressure Drop (bar)")
title("Pressure drop vs diameter ratio")
scatter(dia_rat, delta_P_total)
plot(dia_rat, delta_P_total)
hold on

figure(7)
grid on
xlabel("L_2/D_2")
ylabel("Pressure Drop (bar)")
title("Pressure drop vs length to diameter ratio")
scatter(lengh_dia, delta_P_total)
plot(lengh_dia, delta_P_total)
hold on

figure(8)
grid on
xlabel("Variable (m)")
ylabel("Pressure Drop (bar)")
title("Pressure drop vs Variable")
scatter(variable, delta_P_total)
plot(variable, delta_P_total)
hold on

elseif vary==2 || vary==4 % that is for pipe 2
figure(1)
grid on
xlabel("Variable (m)")
ylabel("Pressure Drop (bar)")
title("Pressure Drop vs Radius")
plot(variable,delta_P12(:,2))
scatter(variable,delta_P12(:,2))
hold on


figure(2)
grid on
xlabel("Variable (m)")
ylabel("Friction Factor")
title("Friction Factor vs Variable")
plot(variable,f12(:,2))
scatter(variable,f12(:,2))
hold on


figure(3)
grid on
xlabel("Reynolds Number")
ylabel("Friction Factor")
title("Reynolds Number vs Friction Factor")
plot(Re12(:,2),f12(:,2))
scatter(Re12(:,2),f12(:,2))
hold on


figure(4)
grid on
xlabel("Variable (m)")
ylabel("Relative Roughness (m^{-1})")
title("Relative Roughness vs Variable")
plot(variable,r_rough12(:,2))
scatter(variable,r_rough12(:,2))
hold on


figure(5)
grid on
xlabel("Variable (m)")
ylabel("Reynolds Number")
title("Reynolds Number vs Variable")
plot(variable,Re12(:,2))
scatter(variable,Re12(:,2))
hold on

figure(6)
grid on
xlabel("D_1/D_2")
ylabel("Pressure Drop (bar)")
title("Pressure drop vs diameter ratio")
scatter(dia_rat, delta_P_total)
plot(dia_rat, delta_P_total)
hold on

figure(7)
grid on
xlabel("L_2/D_2")
ylabel("Pressure Drop (bar)")
title("Pressure drop vs length to diameter ratio")
scatter(lengh_dia, delta_P_total)
plot(lengh_dia, delta_P_total)
hold on

figure(8)
grid on
xlabel("Variable (m)")
ylabel("Pressure Drop (bar)")
title("Pressure drop vs Variable")
scatter(variable, delta_P_total)
plot(variable, delta_P_total)
hold on
end
end
end