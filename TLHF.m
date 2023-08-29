clc;clear all;
close all;
%% input and allotment
set(groot,'defaultLineLineWidth',2.0)

vary = input("What do you want to vary?\nFUEL Parameters:\nA. Fuel Orifice Parameters\n Orifice Radius: 1 \n Number of orifice: 2 \nB. Cross section details: \n channel width: 3 " + ...
    "\n Gap between two channels: 4 \n thetha: 5 \nC. Helical Channel Zone: \n channel height: 6 \n Diameter including channel height (D_sw): 7 \n lead og helix: " + ...
    "8 \n Number of Helical channel: 9 \n Exit diameter: 10 \nOXIDIZER Parameters: \n Offset (Rbx): 11 \n Orifice radius (rbx): 12 \n Element exit radius (rc): 13 \n Number of orifice: 14 \nEnter Your Choice:");
start=input("Enter the starting value of parameter you want to vary in metres: ");
final=input("Enter the ending value of parameter you want to vary in metres: ");
data_points=input("Enter number of data points you want between them: ");
%% independent parameters

% FUEL 
% material
a_rough = linspace(0.03e-3,0.03e-3,data_points); % absolute roughness

% orifice
orifice_rad=linspace(0.9e-3,0.9e-3,data_points);
thickness=linspace(0.9e-3,0.9e-3,data_points);
no_of_orf=linspace(4,4,data_points);

% no channel
no_ch_len=linspace(10.5e-3,10.5e-3,data_points);

%cross section (straight channel)
d1=linspace(1.3e-3,1.3e-3,data_points); % cw: channel width
d2=linspace(1.6e-3,1.6e-3,data_points); % ch: channel height
d3=linspace(1.2e-3,1.2e-3,data_points);
thetha=linspace(15,15,data_points);
no_of_st_ch=linspace(12,12,data_points);
st_ch_len=linspace(10e-3,10e-3,data_points);

%helix
ch=linspace(1.6e-3,1.6e-3,data_points); % ch: channel height
h2=linspace(2.25e-3,2.25e-3,data_points);
r1=linspace(6e-3,6e-3,data_points);
dsw=linspace(11.75e-3,11.75e-3,data_points);
dexit=linspace(10.5e-3,10.5e-3,data_points);
lead=linspace(36e-3,36e-3,data_points);
no_of_hel_ch=linspace(6,6,data_points);
hel_ch_len=linspace(7e-3,7e-3,data_points);

% OXIDIZER
% lox parameters
Rbx=linspace(4.3e-3,4.3e-3,data_points);
rbx=linspace(1.1e-3,1.1e-3,data_points);
rc=linspace(5.5e-3,5.5e-3,data_points); % fixed
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
    d3=linspace(start,final,data_points);
    variable=d3;
elseif vary==5
    thetha=linspace(start,final,data_points);
    variable=thetha;
elseif vary==6
    ch=linspace(start,final,data_points);
    variable=ch;
elseif vary==7
    dsw=linspace(start,final,data_points);
    variable=dsw;
elseif vary==8
    lead=linspace(start,final,data_points);
    variable=lead;
elseif vary==9
    no_of_hel_ch=linspace(start,final,data_points);
    variable=no_of_hel_ch;
elseif vary==10
    dexit=linspace(start,final,data_points);
    variable=dexit;
elseif vary==11
    Rbx=linspace(start,final,data_points);
    variable=Rbx;
elseif vary==12
    rbx=linspace(start, final, data_points);
    variable=rbx;
elseif vary==13
    rc=linspace(start, final, data_points);
    variable=rc;
elseif vary==14
    n=linspace(start, final, data_points);
    variable=n;
end

%% dependent parameters

%definitions
r_rough_orf=a_rough.*0.5./orifice_rad; % relative roughness
d_eff=dsw-ch;
beta=atand(lead./(3.14.*d_eff)); % helix angle
hel_ch_len_eq=hel_ch_len./cosd(beta); % equivalent length

% fuel
mdot_ch4=[0.1072,0.1939,0.223]; %ch4
rho_ch4=linspace(180,180,data_points); %ch4
rho_lox=linspace(1000,1000,data_points);%h20
mu_lox=linspace(853e-6,853e-6,data_points); %h20
mdot_lox=mdot_ch4.*sqrt(rho_lox/rho_ch4); %fuel equivalent water mass flow rate
Temp=linspace(252,252,data_points);

% oxidizer
mdoto=[0.4044,0.7317,0.841];
rho=linspace(1168,1168,data_points);
sigma=linspace(73e-3,73e-3,data_points); % surface tension
%% Calculations
% oxidizer
Ao = Rbx.*rc./(rbx.*rbx.*n); % oxidizer swirl number
orf_area_o=3.14.*rbx.*rbx;
cd_orf=0.69;

for j=1:length(Ao)
fun = @(fio_ox)func(fio_ox,Ao(j)); % calling fsolve function to solve for flow coef
fi0=0.7; % initial fi guess
fio_ox=fsolve(fun,fi0);
fo_ox(j)=fio_ox;
end
sol_to = param_calc_o(fo_ox, Ao, orf_area_o, rho, mdoto,variable,cd_orf,n,rc,data_points);

% fuel

area_cross_sec=0.5.*(d1+d1+2.*ch.*tand(thetha)).*ch;
perimeter_cross_sec=d1+d1+2.*ch.*tand(thetha)+2.*(ch./cosd(thetha));
area_hel_passage=no_of_hel_ch.*area_cross_sec;
A=3.14.*d_eff.*dexit.*cosd(beta)./(4.*area_hel_passage); % Swirl Number

% deltaP calculation due to orifice

velf=zeros(data_points);
for i=1:length(mdot_lox)
mdot_orf=mdot_lox(i)./no_of_orf; % mass flow rate through eachorifice
area_orf=3.14.*orifice_rad.*orifice_rad;
velocity_orf=mdot_orf./(rho_lox.*area_orf);
Re_orf=rho_lox.*velocity_orf.*2.*orifice_rad./mu_lox;
f_orf=friction_factor(Re_orf,r_rough_orf);
deltaP_orf=(mdot_orf./(cd_orf.*area_orf.*sqrt(2.*9.81.*10000.*rho_lox))).^2;
delP_orf_arr(i,:)=deltaP_orf;

% deltaP calculation due to straight channel: PLEASE IGNORE FROM HERE TILL LINE 165
mdot_st_ch=mdot_lox(i)./no_of_st_ch;
area_st_ch=area_cross_sec; 
perimeter_st_ch=perimeter_cross_sec;
D_st_ch=4.*area_st_ch./perimeter_st_ch;
vel_st_ch=mdot_st_ch./(rho_lox.*area_st_ch);
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

% deltaP calculation due to  helical channel
mdot_hel_ch=mdot_lox(i)./no_of_hel_ch;
area_hel_ch=area_cross_sec;
perimeter_hel_ch=perimeter_cross_sec;
D_hel_ch=4.*area_hel_ch./perimeter_hel_ch; % hydraulic diameter
vel_hel_ch=mdot_hel_ch./(rho_lox.*area_hel_ch);
Re_hel_ch=rho_lox.*vel_hel_ch.*D_hel_ch./(mu_lox);
r_rough_hel_ch=a_rough./D_hel_ch;
for j=1:length(r_rough_hel_ch)
fric_hel_ch=@(ff)friction_factor_finder(ff,Re_hel_ch,r_rough_hel_ch);
f0=0.002;
ff=fsolve(fric_hel_ch,f0);
f_hel_ch(j)=ff;
end
deltaP_hel_ch=(f_hel_ch.*hel_ch_len_eq.*vel_hel_ch.*vel_hel_ch./(2.*9.81.*D_hel_ch)).*9.81.*10^(-5).*rho_lox; % pressure drop calculated from flv formula
cd_assumption=0.7;
for j=1:length(A)
fun = @(fioo)func(fioo,A(j)); % calling fsolve function to solve for fi
fi0=0.2; % initial fi guess
fioo=fsolve(fun,fi0);
fio(j)=fioo;
end
sol_tf = param_calc(fio, A, area_hel_ch, rho_lox, mdot_hel_ch,cd_assumption,dexit,no_of_hel_ch,beta);
velf(i,:)=sol_tf(6,:);
deltaP_hel_arr(i,:)=sol_tf(5,:);
deltaP_overall_HF(i,:)=deltaP_hel_arr(i,:)+delP_orf_arr(i,:);


%% plots

    

    figure(16)
    grid on
    title("Fuel Orifice")
    xlabel("Variable");
    ylabel("Del P orifice(bar)");
    scatter(variable,deltaP_orf)
    plot(variable,deltaP_orf)
    hold on

    figure(15)
    grid on
    title("Fuel")
    xlabel("Variable");
    ylabel("Del P overall (bar)");
    scatter(variable,deltaP_overall_HF)
    plot(variable,deltaP_overall_HF)
    hold on

 


    figure(14)
    grid on
    title("Fuel Helical")
    xlabel("Variable");
    ylabel("A");
    scatter(variable,A)
    plot(variable,A)
    hold on

    figure(13)
    grid on
    title("Fuel Helical")
    xlabel("Swirl Number");
    ylabel("Del P (helical)");
    scatter(A,sol_tf(3,:))
    plot(A,sol_tf(3,:))
    hold on

    figure(12)
    grid on
    title("Fuel Helical")
    xlabel("Swirl Number");
    ylabel("Del P (Swirl)");
    scatter(A,sol_tf(4,:))
    plot(A,sol_tf(4,:))
    hold on

    figure(11)
    grid on
    title("Fuel Helical")
    xlabel("Swirl Number");
    ylabel("Del P (swirl plus helical)");
    scatter(A,sol_tf(5,:))
    plot(A,sol_tf(5,:))
    hold on

   

    figure(10)
    grid on
    title("Fuel Helical")
    xlabel("A");
    ylabel("cd,fi,full cone angle (*100)");
    scatter(A,sol_tf(1,:))
    plot(A,sol_tf(1,:))
    scatter(A,sol_tf(2,:)*0.01)
    plot(A,sol_tf(2,:)*0.01)
    scatter(A,fio)
    plot(A,fio)
    hold on


    figure(9)
    grid on
    title("Fuel Helical")
    xlabel("A");
    ylabel("Velocity-helical channel axial");
    scatter(A,sol_tf(6,:))
    plot(A,sol_tf(6,:))
    hold on 

    figure(91)
    grid on
    title("Fuel Helical")
    xlabel("A");
    ylabel("Velocity-exit axial");
    scatter(A,sol_tf(7,:))
    plot(A,sol_tf(7,:))
    hold on 

    figure(92)
    grid on
    title("Fuel Helical")
    xlabel("A");
    ylabel("resultant velocity in helical channel");
    scatter(A,sol_tf(8,:))
    plot(A,sol_tf(8,:))
    hold on 

    figure(93)
    grid on
    title("Fuel Helical")
    xlabel("A");
    ylabel("element exit resultant velocity");
    scatter(A,sol_tf(9,:))
    plot(A,sol_tf(9,:))
    hold on 
end


for q=1:length(mdoto)
    vel_ox(q,:)=sol_to(q+6,:);
    VR(q,:)=velf(q,:)./vel_ox(q,:); % velocity ratio
    MR(q,:)=(rho_lox.*velf(q,:).*velf(q,:))./(rho.*vel_ox(q,:).*vel_ox(q,:)); % momentum ratio
     figure(8)
    grid on
    title("Velocity Ratio")
    xlabel("variable");
    ylabel("VR");
    scatter(variable,VR(q,:))
    plot(variable,VR(q,:))
    hold on 

    figure(7)
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

function calc=param_calc(fi, A, area_hel_ch, rho_lox, mdot_hel_ch,cd_assumption,de,n,beta) % for fuel
cd=sqrt(fi.^3./(2-fi));
alpha = 2.*atand((2.*cd.*A)./(sqrt((1+sqrt(1-fi)).^2-4.*cd.*cd.*A.*A)));
deltaP_hel_ch_1=(mdot_hel_ch./(cd_assumption.*area_hel_ch.*sqrt(2.*9.81.*rho_lox.*10000))).^2;
deltaP_hel_ch_2=((n.*mdot_hel_ch)./(cd.*3.14.*de.*de.*0.25.*sqrt(2.*9.81.*rho_lox.*10000))).^2;
total_delp_helical=deltaP_hel_ch_2+deltaP_hel_ch_1;
hel_lox_axial_vel=(mdot_hel_ch./(rho_lox.*area_hel_ch)).*cosd(beta); % hel ch axial vel
hel_lox_tang=hel_lox_axial_vel.*tand(alpha./2); % tangential velocity
hel_lox_res=sqrt(hel_lox_tang.*hel_lox_tang + hel_lox_axial_vel.*hel_lox_axial_vel);
exit_lox_axial_vel=((mdot_hel_ch.*n)./(rho_lox.*0.25.*3.14.*de.*de.*fi)); % exit axial vel
exit_lox_tang=exit_lox_axial_vel.*tand(alpha./2);
exit_lox_res=sqrt(exit_lox_tang.*exit_lox_tang + exit_lox_axial_vel.*exit_lox_axial_vel);
calc=[cd;alpha; deltaP_hel_ch_1;deltaP_hel_ch_2;total_delp_helical; hel_lox_axial_vel; exit_lox_axial_vel;hel_lox_res;exit_lox_res];
end

function calc=param_calc_o(fi, A, orifice_area, rho_lox, mdot_lox,variable,cd_orf,n,rc,data_points) % for oxidizer
velo=zeros(data_points);
delp_tot_ox=zeros(data_points);
cda=zeros(data_points);
for i=1:length(mdot_lox)
cd=sqrt(fi.^3./(2-fi));
alpha = 2.*atand((2.*cd.*A)./(sqrt((1+sqrt(1-fi)).^2-4.*cd.*cd.*A.*A)));
delp1=(mdot_lox(i)./(cd_orf.*n.*orifice_area.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to orifice
delp2=(mdot_lox(i)./(cd.*3.14.*rc.*rc.*sqrt(2.*9.81.*rho_lox.*10000))).^2; % due to swirl with effects of friction included in cd
delp_tot=delp2+delp1;
lox_axial_vel=mdot_lox(i)./(rho_lox.*3.14.*rc.*rc.*fi);
cda(i,:)=mdot_lox(i)./sqrt(2.*rho_lox.*delp_tot.*10^5);
velo(i,:)=lox_axial_vel;
delp_tot_ox(i,:)=delp_tot;
calc=[cd;alpha;delp1;delp2;delp_tot;lox_axial_vel;velo(1,:);velo(2,:);velo(3,:);delp_tot_ox(1,:);delp_tot_ox(2,:);delp_tot_ox(3,:);cda(1,:)];


figure(6)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("deltaP1 (due to orifice)")
plot(A,delp1)
scatter(A,delp1)
hold on

figure(5)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("deltaP2 (due to swirl)")
plot(A,delp2)
scatter(A,delp2)
hold on

figure(4)
grid on
title("Oxidizer")
ylabel("Swirl Number")
xlabel("Variable")
plot(variable, A)
scatter(variable, A)
hold on

figure(3)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("deltaP_total (orifice plus swirl)")
plot(A,delp_tot)
scatter(A,delp_tot)
hold on

figure(2)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("velocity (m/s)")
plot(A,lox_axial_vel)
scatter(A,lox_axial_vel)
hold on

figure(1)
grid on
title("Oxidizer")
xlabel("Swirl Number")
ylabel("cd,fi,full cone angle(*100")
plot(A,fi)
scatter(A,fi)
plot(A,cd)
scatter(A,cd)
plot(A,alpha*0.01)
scatter(A,alpha*0.01)
hold on
end
end