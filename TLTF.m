% lox swirl analysis using tangential holes
clc;clear all;close all;

set(groot,'defaultLineLineWidth',2.0)

vary = input("What do you want to vary? \n Fuel Parameters:\n 1. Offset: \n 2. Orifice radius: \n 3. Element exit radius: " + ...
    "\n 4. Number of orifice: \n 5. Gap 1: \n 6. Gap 2: \n Oxidizer Parameter: \n 7. Element exit radius: \n 8. Orifice Radius: \n 9. Offset: " + ...
    "\n 10. Number of orifice: \n 11. Lip: \n Enter your choice: ");
start=input("Enter the starting value of parameter you want to vary in metres: ");
final=input("Enter the ending value of parameter you want to vary in metres: ");
data_points=input("Enter number of data points you want between them: ");

% overall parameters
% oxidizer
d1=linspace(5.5e-3,5.5e-3,data_points);
d2=linspace(3e-3,3e-3,data_points);
l1=linspace(14e-3,14e-3,data_points);
l2=linspace(36.2e-3,36.2e-3,data_points);

% fuel
lip=linspace(0.9e-3,0.9e-3,data_points);
e=linspace(3.4e-3,3.4e-3,data_points);
a_gap2=linspace(0.4e-3,0.4e-3,data_points);
a_gap1=linspace(0.8e-3,0.8e-3,data_points);
lf1=linspace(20.8e-3,20.8e-3,data_points);
lf2=linspace(9.7e-3,9.7e-3,data_points);


% ch4 parameters water equivalent
Rbxf=linspace(3.4e-3,3.4e-3,data_points); % offset
rbxf=linspace(2e-3,2e-3,data_points); % orf dia
nf=linspace(4,4,data_points);
mdot_ch4=[0.1072, 0.1939, 0.223];
mdot_f=mdot_ch4.*sqrt(1000/180);
rho_f=linspace(1000,1000,data_points);
mu_f=linspace(853e-6,853e-6,data_points);
sigma_f=linspace(853e-6,853e-6,data_points);

% lox parameters
Rbxo=linspace(3.4e-3,3.4e-3,data_points); % offset
rbxo=linspace(2e-3,2e-3,data_points); % orf dia
rco=linspace(5.5e-3,5.5e-3,data_points); % element exit dia
no=linspace(4,4,data_points);
mdot_o=[0.4044,0.7317,0.841];
rho_o=linspace(1168,1168,data_points);
mu_o=linspace(215.39e-6,215.39e-6,data_points);
sigma_o=linspace(73e-3,73e-3,data_points);

if vary==1
    Rbxf=linspace(start,final,data_points);
    variable=Rbxf;
elseif vary==2
    rbxf=linspace(start, final, data_points);
    variable=rbxf;
elseif vary==4
    nf=linspace(start, final, data_points);
    variable=nf; 
elseif vary==5
    a_gap1=linspace(start, final, data_points);
    variable=a_gap1;
elseif vary==6
    a_gap2=linspace(start, final, data_points);
    variable=a_gap2;
elseif vary==7
    rco=linspace(start, final, data_points);
    variable=rco;
elseif vary==8
    rbxo=linspace(start, final, data_points);
    variable=rbxo; 
elseif vary==9
    Rbxo=linspace(start, final, data_points);
    variable=Rbxo;
elseif vary==10
    no=linspace(start, final, data_points);
    variable=no;
end

df1=d1+lip; % inner radius
df2=d2+lip;
df11=df2+a_gap1; % outer radius
df22=df2+a_gap2;
rcf=df11; % element exit dia

% f
exit_area_f=3.14.*(df11.*df11-df2.*df2);
exit_perimeter_f=2.*3.14.*(df11+df2);
hydraulic_f=4.*exit_area_f./exit_perimeter_f;
Af = Rbxf.*rcf./(rbxf.*rbxf.*nf);
orf_area_f=3.14.*rbxf.*rbxf;
cd_orf_f=0.75;

% o
exit_area_o=3.14.*rco.*rco;
exit_perimeter_o=2.*3.14.*rco;
hydraulic_o=4.*exit_area_o./exit_perimeter_o;
Ao = Rbxo.*rco./(rbxo.*rbxo.*no);
orf_area_o=3.14.*rbxo.*rbxo;
cd_orf_o=0.75;

% f
for j=1:length(Af)
funf = @(fio)func(fio,Af(j)); % calling fsolve function to solve for fi
fi0=0.7; % initial fi guess
fio=fsolve(funf,fi0);
ff(j)=fio;
end
sol_tf = param_calc_f(ff, Af, orf_area_f, rho_f, mdot_f,variable,cd_orf_f,nf,data_points,df2,df11);

% o
for j=1:length(Ao)
funo = @(fio)func(fio,Ao(j)); % calling fsolve function to solve for fi
fi0=0.7; % initial fi guess
fio=fsolve(funo,fi0);
fo(j)=fio;
end
sol_to = param_calc_o(fo, Ao, orf_area_o, rho_o, mdot_o,variable,cd_orf_o,no,data_points,rco);

% Velocity Ratio
for k=1:length(mdot_o)
figure(17)
grid on
xlabel("Variable")
ylabel("velocity Ratio")
plot(variable,sol_tf(7,:)./sol_to(7,:))
scatter(variable,sol_tf(7,:)./sol_to(7,:))
plot(variable,sol_tf(8,:)./sol_to(8,:))
scatter(variable,sol_tf(8,:)./sol_to(8,:))
plot(variable,sol_tf(9,:)./sol_to(9,:))
scatter(variable,sol_tf(9,:)./sol_to(9,:))
hold on

figure(18)
grid on
xlabel("Variable")
ylabel("Momentum Ratio")
plot(variable,(rho_f.*sol_tf(7,:).*sol_tf(7,:))./(rho_o.*sol_to(7,:).*sol_to(7,:)))
scatter(variable,(rho_f.*sol_tf(7,:).*sol_tf(7,:))./(rho_o.*sol_to(7,:).*sol_to(7,:)))
plot(variable,(rho_f.*sol_tf(8,:).*sol_tf(8,:))./(rho_o.*sol_to(8,:).*sol_to(8,:)))
scatter(variable,(rho_f.*sol_tf(8,:).*sol_tf(8,:))./(rho_o.*sol_to(8,:).*sol_to(8,:)))
plot(variable,(rho_f.*sol_tf(9,:).*sol_tf(9,:))./(rho_o.*sol_to(9,:).*sol_to(9,:)))
scatter(variable,(rho_f.*sol_tf(9,:).*sol_tf(9,:))./(rho_o.*sol_to(9,:).*sol_to(9,:)))
hold on

figure(19)
grid on
xlabel("Variable")
ylabel("Oxidizer Reynolds Number")
plot(variable,sol_to(7,:).*rho_o.*hydraulic_o./mu_o);
scatter(variable,sol_to(7,:).*rho_o.*hydraulic_o./mu_o)
plot(variable,sol_to(8,:).*rho_o.*hydraulic_o./mu_o)
scatter(variable,sol_to(8,:).*rho_o.*hydraulic_o./mu_o)
plot(variable,sol_to(9,:).*rho_o.*hydraulic_o./mu_o)
scatter(variable,sol_to(9,:).*rho_o.*hydraulic_o./mu_o)
hold on

figure(20)
grid on
xlabel("Variable")
ylabel("Fuel Reynolds Number")
plot(variable,sol_tf(7,:).*rho_f.*hydraulic_f./mu_f);
scatter(variable,sol_tf(7,:).*rho_f.*hydraulic_f./mu_f)
plot(variable,sol_tf(8,:).*rho_f.*hydraulic_f./mu_f)
scatter(variable,sol_tf(8,:).*rho_f.*hydraulic_f./mu_f)
plot(variable,sol_tf(9,:).*rho_f.*hydraulic_f./mu_f)
scatter(variable,sol_tf(9,:).*rho_f.*hydraulic_f./mu_f)
hold on

figure(21)
grid on
xlabel("Variable")
ylabel("Fuel Weber Number")
plot(variable,sol_tf(7,:).*sol_tf(7,:).*rho_f.*hydraulic_f./sigma_f);
scatter(variable,sol_tf(7,:).*sol_tf(7,:).*rho_f.*hydraulic_f./sigma_f)
plot(variable,sol_tf(8,:).*sol_tf(8,:).*rho_f.*hydraulic_f./sigma_f)
scatter(variable,sol_tf(8,:).*sol_tf(8,:).*rho_f.*hydraulic_f./sigma_f)
plot(variable,sol_tf(9,:).*sol_tf(9,:).*rho_f.*hydraulic_f./sigma_f)
scatter(variable,sol_tf(9,:).*sol_tf(9,:).*rho_f.*hydraulic_f./sigma_f)
hold on

figure(22)
grid on
xlabel("Variable")
ylabel("Oxidizer Weber Number")
plot(variable,sol_to(7,:).*sol_to(7,:).*rho_o.*hydraulic_o./sigma_o);
scatter(variable,sol_to(7,:).*sol_to(7,:).*rho_o.*hydraulic_o./sigma_o)
plot(variable,sol_to(8,:).*sol_to(8,:).*rho_o.*hydraulic_o./sigma_o)
scatter(variable,sol_to(8,:).*sol_to(8,:).*rho_o.*hydraulic_o./sigma_o)
plot(variable,sol_to(9,:).*sol_to(9,:).*rho_o.*hydraulic_o./sigma_o)
scatter(variable,sol_to(9,:).*sol_to(9,:).*rho_o.*hydraulic_o./sigma_o)
hold on
end


% functions
function F = func(f,A)
F=1-sqrt((A.^2./(1-f))+(f.^(-2))).*sqrt(f.^3/(2-f));
end

function calc=param_calc_f(fi, A, orifice_area, rho_lox, mdot_lox,variable,cd_orf,nf,data_points,df2,df11)
vel=zeros(data_points);
for i=1:length(mdot_lox)
cd=sqrt(fi.^3./(2-fi));
alpha = 2.*atand((2.*cd.*A)./(sqrt((1+sqrt(1-fi)).^2-4.*cd.*cd.*A.*A)));
delp1=(mdot_lox(i)./(cd_orf.*nf.*orifice_area.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to orifice
delp2=(mdot_lox(i)./(cd.*3.14.*(df11.*df11-df2.*df2).*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to swirl with effects of friction included in cd
delp_tot=delp2+delp1;
lox_axial_vel=mdot_lox(i)./(rho_lox.*3.14.*(df11.*df11-df2.*df2).*fi);
vel(i,:)=lox_axial_vel;
calc=[cd;alpha; delp1; delp2; delp_tot;lox_axial_vel;vel(1,:);vel(2,:);vel(3,:)];



figure(3)
grid on
title("FUEL")
xlabel("Swirl Number")
ylabel("fi,cd,angle")
plot(A,fi)
scatter(A,fi)
plot(A,alpha.*0.01)
scatter(A,alpha.*0.01)
plot(A,cd)
scatter(A,cd)
hold on

figure(4)
grid on
title("FUEL")
xlabel("Swirl Number")
ylabel("deltaP1 (orf)")
plot(A,delp1)
scatter(A,delp1)
hold on

figure(5)
grid on
title("FUEL")
xlabel("Swirl Number")
ylabel("deltaP2 (Swirl)")
plot(A,delp2)
scatter(A,delp2)
hold on

figure(6)
grid on
title("FUEL")
ylabel("Swirl Number")
xlabel("Variable")
plot(variable, A)
scatter(variable, A)
hold on

figure(7)
grid on
title("FUEL")
xlabel("Swirl Number")
ylabel("deltaP_total")
plot(A,delp_tot)
scatter(A,delp_tot)
hold on

figure(8)
grid on
title("FUEL")
xlabel("Swirl Number")
ylabel("velocity (m/s)")
plot(A,lox_axial_vel)
scatter(A,lox_axial_vel)
hold on
end
end

function calc=param_calc_o(fi, A, orifice_area, rho_lox, mdot_lox,variable,cd_orf,n,data_points,rc)
velo=zeros(data_points);
for i=1:length(mdot_lox)
cd=sqrt(fi.^3./(2-fi));
alpha = 2.*atand((2.*cd.*A)./(sqrt((1+sqrt(1-fi)).^2-4.*cd.*cd.*A.*A)));
delp1=(mdot_lox(i)./(cd_orf.*n.*orifice_area.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to orifice
delp2=(mdot_lox(i)./(cd.*3.14.*rc.*rc.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to swirl with effects of friction included in cd
delp_tot=delp2+delp1;
lox_axial_vel=mdot_lox(i)./(rho_lox.*3.14.*rc.*rc.*fi);
velo(i,:)=lox_axial_vel;
calc=[cd;alpha; delp1;delp2; delp_tot;lox_axial_vel;velo(1,:);velo(2,:);velo(3,:)];



figure(10)
grid on
xlabel("Swirl Number")
ylabel("cd,fi,Cone angle")
plot(A,alpha)
scatter(A,alpha)
plot(A,fi)
scatter(A,fi)
plot(A,cd)
scatter(A,cd)
hold on


figure(12)
grid on
xlabel("Swirl Number")
ylabel("deltaP1 (orf)")
plot(A,delp1)
scatter(A,delp1)
hold on

figure(13)
grid on
xlabel("Swirl Number")
ylabel("deltaP2 (swirl)")
plot(A,delp2)
scatter(A,delp2)
hold on

figure(14)
grid on
ylabel("Swirl Number")
xlabel("Variable")
plot(variable, A)
scatter(variable, A)
hold on

figure(15)
grid on
xlabel("Swirl Number")
ylabel("deltaP_total")
plot(A,delp_tot)
scatter(A,delp_tot)
hold on

figure(16)
grid on
xlabel("Swirl Number")
ylabel("velocity (m/s)")
plot(A,lox_axial_vel)
scatter(A,lox_axial_vel)
hold on
end
end