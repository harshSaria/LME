clc;clear all;
close all;
%% input and allotment
set(groot,'defaultLineLineWidth',2.0)

vary = input("What do you want to vary?\nFUEL Parameters:\n A. Fuel Orifice Parameters\n Orifice Radius: 1 \nOrifice length: 2 \nNumber of orifice: 3 \n B. No channel Zone\n Length of no channel zone: 4 \nC. Straight channel Zone: \nd1: 5 " + ...
    "\nd2: 6 \nd3: 7 \nthetha: 8 \nNumber of staright channels: 9 \nLength of staright channels: 10 \nOXIDIZER Parameters: \nRbx: 11 \nrbx: 12 \nrc: 13 \nn: 14 \nEnter Your Choice:");
start=input("Enter the starting value of parameter you want to vary in metres: ");
final=input("Enter the ending value of parameter you want to vary in metres: ");
data_points=input("Enter number of data points you want between them: ");
%% independent parameters

% FUEL 
% material
a_rough = linspace(0.03e-3,0.03e-3,data_points);

% orifice
orifice_rad=linspace(1e-3,1e-3,data_points);
thickness=linspace(0.9e-3,0.9e-3,data_points);
no_of_orf=linspace(4,4,data_points);

% no channel
no_ch_len=linspace(10.5e-3,10.5e-3,data_points);

%cross section (straight channel)
d1=linspace(1.8e-3,1.8e-3,data_points); %cw
d2=linspace(1.2e-3,1.2e-3,data_points); %ch
d3=linspace(1.2e-3,1.2e-3,data_points);
thetha=linspace(15,15,data_points);
no_of_st_ch=linspace(6,6,data_points);
st_ch_len=linspace(35e-3,35e-3,data_points);

% OXIDIZER
% lox parameters
Rbx=linspace(3.4e-3,3.4e-3,data_points);
rbx=linspace(1.2e-3,1.2e-3,data_points);
rc=linspace(5.5e-3,5.5e-3,data_points);
n=linspace(4,4,data_points);

%% re-allotment

if vary==1
    orifice_rad=linspace(start,final,data_points);
    variable=orifice_rad;
elseif vary==2
    no_of_orf=linspace(start,final,data_points);
    variable=no_of_orf;
elseif vary==3
    d1=linspace(start,final,data_points);
    variable=d1;
elseif vary==4
    d2=linspace(start,final,data_points);
    variable=d2;
elseif vary==5
    d3=linspace(start,final,data_points);
    variable=d3;
elseif vary==6
    thetha=linspace(start,final,data_points);
    variable=thetha;
elseif vary==7
    no_of_st_ch=linspace(start,final,data_points);
    variable=no_of_st_ch;
elseif vary==8
    st_ch_len=linspace(start,final,data_points);
    variable=st_ch_len;
elseif vary==9
    Rbx=linspace(start,final,data_points);
    variable=Rbx;
elseif vary==10
    rbx=linspace(start, final, data_points);
    variable=rbx;
elseif vary==11
    rc=linspace(start, final, data_points);
    variable=rc;
elseif vary==12
    n=linspace(start, final, data_points);
    variable=n;
end

%% dependent parameters

%definitions
r_rough_orf=a_rough.*0.5./orifice_rad;

% fuel
mdot_ch4=[0.1072,0.1939,0.223];
rho_ch4=linspace(180,180,data_points);
rho_lox=linspace(1000,1000,data_points);
mdot_lox=mdot_ch4.*sqrt(rho_lox/rho_ch4);
mu_lox=linspace(853e-6,853e-6,data_points);
sigma_lox=linspace(1,1,data_points);
Temp=linspace(252,252,data_points);

% oxidizer
mdoto=[0.4044,0.7317,0.841];
rho=linspace(1168,1168,data_points);
sigma=linspace(73e-3,73e-3,data_points);

%% Calculations
% oxidizer
Ao = Rbx.*rc./(rbx.*rbx.*n);
orf_area_o=3.14.*rbx.*rbx;
cd_orf=0.75;

for j=1:length(Ao)
fun = @(fio_ox)func(fio_ox,Ao(j)); % calling fsolve function to solve for fi
fi0=0.7; % initial fi guess
fio_ox=fsolve(fun,fi0);
fo_ox(j)=fio_ox;
end
sol_to = param_calc_o(fo_ox, Ao, orf_area_o, rho, mdoto,variable,cd_orf,n,rc,data_points);

% fuel
% cross section
area_cross_sec=0.5.*(d1+d1+2.*d2.*tand(thetha)).*d2;
perimeter_cross_sec=d1+d1+2.*d2.*tand(thetha)+2.*(d2./cosd(thetha));

% deltaP calculation due to orifice
k=0.1;
velf=zeros(data_points);
for i=1:length(mdot_lox)
mdot_orf=mdot_lox(i)./no_of_orf;
area_orf=3.14.*orifice_rad.*orifice_rad;
velocity_orf=mdot_orf./(rho_lox.*area_orf);
Re_orf=rho_lox.*velocity_orf.*2.*orifice_rad./mu_lox;
f_orf=friction_factor(Re_orf,r_rough_orf);
deltaP_orf=(mdot_orf./(0.75.*area_orf.*sqrt(2.*9.81.*10000.*rho_lox))).^2;
delP_orf_arr(i,:)=deltaP_orf;

% deltaP calculation due to straight channel
mdot_st_ch=mdot_lox(i)./no_of_st_ch;
area_st_ch=area_cross_sec; 
perimeter_st_ch=perimeter_cross_sec;
D_st_ch=4.*area_st_ch./perimeter_st_ch;
vel_st_ch=mdot_st_ch./(rho_lox.*area_st_ch);
velf(i,:)=vel_st_ch;
Re_st_ch=rho_lox.*vel_st_ch.*D_st_ch./(mu_lox);
r_rough_st_ch=a_rough./D_st_ch;
for j=1:length(r_rough_st_ch)
fric=@(f)friction_factor_finder(f,Re_st_ch,r_rough_st_ch);
f0=0.002;
f=fsolve(fric,f0);
f_st_ch(j)=f;
end
deltaP_st_ch=(f_st_ch.*st_ch_len.*vel_st_ch.*vel_st_ch./(2.*9.81.*D_st_ch)).*rho_lox.*9.81.*10^(-5);
delp_st_ch_arr(i,:)=deltaP_st_ch;

deltaP_overall_SCF(i,:)=delp_st_ch_arr(i,:)+delP_orf_arr(i,:);


%% plots

  

    figure(2)
    grid on
    title("Fuel Orifice")
    xlabel("Variable");
    ylabel("Del P (bar)");
    scatter(variable,deltaP_orf)
    plot(variable,deltaP_orf)
    hold on

    figure(3)
    grid on
    title("Fuel straight channel")
    xlabel("Variable");
    ylabel("Del P");
    scatter(variable,deltaP_st_ch)
    plot(variable,deltaP_st_ch)
    hold on

      figure(4)
    grid on
    title("Fuel")
    xlabel("Variable");
    ylabel("Del P overall - orf plus st ch(bar)");
    scatter(variable,deltaP_overall_SCF)
    plot(variable,deltaP_overall_SCF)
    hold on
end

for q=1:length(mdoto)
    vel_ox(q,:)=sol_to(q+6,:);
    VR(q,:)=velf(q,:)./vel_ox(q,:);
    MR(q,:)=(rho_lox.*velf(q,:).*velf(q,:))./(rho.*vel_ox(q,:).*vel_ox(q,:));
     figure(5)
    grid on
    title("Velocity Ratio")
    xlabel("variable");
    ylabel("VR");
    scatter(variable,VR(q,:))
    plot(variable,VR(q,:))
    hold on 

    figure(6)
    grid on
    title("Momentum Ratio")
    xlabel("variable");
    ylabel("MR");
    scatter(variable,MR(q,:))
    plot(variable,MR(q,:))
    hold on 
end


%% functions

function f = friction_factor(Re, r_rough)
f=(-2.*log10((r_rough./3.7)-(5.02./Re).*log10(r_rough-(5.02./Re).*((r_rough./3.7)+(13./Re))))).^(-2);
% f=(-1.8.*log10((r_rough./3.7).^1.11 + (6.9./Re))).^(-2);
% f=(1.14-2.*log10(r_rough+(21.25./Re.^0.9))).^(-2);
% f=(-2.*log10((r_rough./3.7)+(5.74./Re.^0.9))).^(-2);
end

function F = friction_factor_finder(f, Re, r_rough)
F=f-(-2.*log10((r_rough./3.7)-(5.02./Re).*log10(r_rough-(5.02./Re).*((r_rough./3.7)+(13./Re))))).^(-2);
% f=(-1.8.*log10((r_rough./3.7).^1.11 + (6.9./Re))).^(-2);
% f=(1.14-2.*log10(r_rough+(21.25./Re.^0.9))).^(-2);
% f=(-2.*log10((r_rough./3.7)+(5.74./Re.^0.9))).^(-2);
end

function F = func(f,A)
F=1-sqrt((A.^2/(1-f))+(f.^(-2))).*sqrt(f.^3./(2-f));
end


function calc=param_calc_o(fi, A, orifice_area, rho_lox, mdot_lox,variable,cd_orf,n,rc,data_points)
velo=zeros(data_points);
delp_tot_ox=zeros(data_points);
for i=1:length(mdot_lox)
cd=sqrt(fi.^3./(2-fi));
alpha = 2.*atand((2.*cd.*A)./(sqrt((1+sqrt(1-fi)).^2-4.*cd.*cd.*A.*A)));
delp1=(mdot_lox(i)./(cd_orf.*n.*orifice_area.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to orifice
delp2=(mdot_lox(i)./(cd.*3.14.*rc.*rc.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to swirl with effects of friction included in cd
delp_tot=delp2+delp1;
lox_axial_vel=mdot_lox(i)./(rho_lox.*3.14.*rc.*rc.*fi);
velo(i,:)=lox_axial_vel;
delp_tot_ox(i,:)=delp_tot;
calc=[cd;alpha;delp1;delp2;delp_tot;lox_axial_vel;velo(1,:);velo(2,:);velo(3,:);delp_tot_ox(1,:);delp_tot_ox(2,:);delp_tot_ox(3,:)];

figure(24)
grid on
title("Oxidizer")
xlabel("A")
ylabel("cd fi alpha(*100)")
plot(A,fi)
scatter(A,fi)
plot(A,alpha.*0.01)
scatter(A,alpha.*0.01)
plot(A,cd)
scatter(A,cd)
hold on

figure(20)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("deltaP1 (due to orifice)")
plot(A,delp1)
scatter(A,delp1)
hold on

figure(21)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("deltaP2 (due to swirl)")
plot(A,delp2)
scatter(A,delp2)
hold on

figure(22)
grid on
title("Oxidizer")
ylabel("Swirl Number")
xlabel("Variable")
plot(variable, A)
scatter(variable, A)
hold on

figure(23)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("deltaP_total (orifice plus swirl)")
plot(A,delp_tot)
scatter(A,delp_tot)
hold on

figure(24)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("velocity (m/s)")
plot(A,lox_axial_vel)
scatter(A,lox_axial_vel)
hold on


end
end