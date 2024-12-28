% see Shen et al, 2024 for details
% sample code for processing strain stress data
% process S plug
% this function could benifit from the usage of symolic math package(e.g,
% use command: digits(50)

clear all
close all

theta = 45*pi/180;
theta_f = 37.1*pi/180
tau1 = 60*pi/180
tau2 = 60*pi/180

eV1_fail = 1 % set flag that ev1 has failed

% from 1st to 5th colum, present epV1 ef45, eh1, eh2 and eh3, see shen et al.
strain_data = [
    1207	1348	1716	2337	181
    1006	1366	1712	2360	249
    1039	1317	1540	2288	-256
    784     1268	1532	2326	-196
    863     1211	1293	2216	-560
    636     1158	1272	2244	-533
    890     1079	938     2035	-996
    281     927     950     2132	-909
    427     888     666     1942	-1286
    142     797     672     2001	-1242
    212     732     438     1807	-1594
    -299	581     503     1950	-1504
    -172	472     106     1592	-2120
    -9999	255     183     1805	-1984
    -9999	109     -384	1305	-2740
];
% NOTE  -9999 means strain gage failes
% first column: axial pressure; second column: radial pressure
pressure_data = [
0	0
2	0
2	2
5	2
5	5
10	5
10	10
15	10
15	15
20	15
20	20
30	20
30	30
45	30
45	45
];


eV1 = strain_data(1:end,1)-strain_data(1,1); % fail
e45 = strain_data(1:end,2)-strain_data(1,2); 
eH1 = strain_data(1:end,3)-strain_data(1,3);
eH2 = strain_data(1:end,4)-strain_data(1,4);
eH3 = strain_data(1:end,5)-strain_data(1,5); 


eV1 = eV1/1e6*-1;
e45 = e45/1e6*-1;
eH1 = eH1/1e6*-1;
eH2 = eH2/1e6*-1;
eH3 = eH3/1e6*-1;

axial_press = pressure_data(1:end, 1);
radial_press = pressure_data(1:end, 2);

figure
subplot(121)
plot(axial_press, e45, axial_press, eH1, axial_press,eH2,axial_press, eV1, axial_press, eH3)

if eV1_fail ==1
    for i = 1:length(eH1)
        %function [ep1, ep3] = V_strain_fix(ep1,ep_theta,ep3, theta)
        [eV1_c(i)] = V_strain_fix(NaN,e45(i),eH1(i), theta); % ev1 fail
    end
else
    eV1_c = eV1;
end



% from 30 -40 MPa, use eV1_c
eV1_c = eV1_c - (eV1_c(end-2) - eV1(end-2))
for i = 1:length(eH1)
    if i < length(eH1)-2
        eV1_c(i)=eV1(i)
    end
end

subplot(122)
plot(axial_press,eV1, axial_press,eV1_c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S plug process
%

for i = 1:length(eH1) 
    [e2p3(i), e3p3(i),om(i)] = H_strain_gages_inversion(eH1(i), eH2(i), eH3(i),NaN,tau1,tau2); % find what are eH1 and eH2
end
figure
plot(axial_press,om/pi*180)
title('Omega angle')



% obtain ep1p, ep2p, ep3p

e2p = e2p3;
e3p = e3p3;
e1p = eV1_c;


%figure
%plot(axial_press,e1p,axial_press,e2p,axial_press,e3p)

%calc S44
sq_alpha = sin(theta_f);
sq_beta = cos(theta_f);
for i = 1:length(e2p)
    S44up(i) = e3p(i)*sq_alpha*sq_beta - e1p(i)*sq_alpha*sq_beta;
    S44down(i) = radial_press(i)*sq_alpha*sq_beta - axial_press(i)*sq_alpha*sq_beta;
    if i > 1
        S44(i) = (S44up(i) - S44up(i-1))/(S44down(i) - S44down(i-1));
    end
    
end



alpha = (sin(theta_f))^2 ; % 
beta = (cos(theta_f))^2;

e2 = e2p;
e1 = e1p(:)*beta + e3p(:)*alpha;
e3 = e1p(:)*alpha + e3p(:)*beta;



R1= [];
R2 = [];
R3 = [];

for i = [3 5 7 9 11 13 15]
    R1(end+1) = (e1(i) - e1(i-1))/(radial_press(i) - radial_press(i-1));
    R2(end+1) = (e2(i) - e2(i-1))/(radial_press(i) - radial_press(i-1));
    R3(end+1) = (e3(i) - e3(i-1))/(radial_press(i) - radial_press(i-1));
end


A1= [];
A2 = [];
A3 = [];

pressure_for_S = []
for i = [2 4 6 8 10 12 14]
    A1(end+1) = (e1(i) - e1(i-1))/(axial_press(i) - axial_press(i-1));
    A2(end+1) = (e2(i) - e2(i-1))/(axial_press(i) - axial_press(i-1));
    A3(end+1) = (e3(i) - e3(i-1))/(axial_press(i) - axial_press(i-1));
    pressure_for_S(end+1)= axial_press(i)
end


a = alpha;
b = beta;
for i = 2:length(A3) %
    
    bxa = [R1(i);R2(i);R3(i);A1(i);A2(i);A3(i)]
    Axa = [a 1 b 0;1 a b 0; 0 0 (1+a) b;b 0 a 0; 0 b a 0; 0 0 b a]
    xout = Axa\bxa;
    S11H(i)= xout(1);
    S12H(i) = xout(2);
    S13H(i) = xout(3);
    S33H(i) = xout(4);
    
    bxa = [R1(i-1);R2(i-1);R3(i-1);A1(i);A2(i);A3(i)]
    AxaA = [a 1 b 0;1 a b 0; 0 0 (1+a) b; b 0 a 0; 0 b a 0; 0 0 b a]
    xout = Axa\bxa;
    S11L(i)= xout(1);
    S12L(i) = xout(2);
    S13L(i) = xout(3);
    S33L(i) = xout(4);
    
    
end
figure
plot(1:length(S11H),S11H,1:length(S11H),S12H,1:length(S11H),S13H,1:length(S11H),S33H)
hold on
plot(1:length(S44),S44)
figure

hold on
pgonS11 = polyshape([pressure_for_S fliplr(pressure_for_S)],[S11H fliplr(S11L)] )
pgonS12 = polyshape([pressure_for_S fliplr(pressure_for_S)],[S12H fliplr(S12L)] )
pgonS13 = polyshape([pressure_for_S fliplr(pressure_for_S)],[S13H fliplr(S13L)] )
pgonS33 = polyshape([pressure_for_S fliplr(pressure_for_S)],[S33H fliplr(S33L)] )
plot(pgonS11)
plot(pgonS12)
plot(pgonS13)
plot(pgonS33)
plot(pressure_for_S,(S11H+S11L)/2,pressure_for_S,(S12H+S12L)/2,pressure_for_S,(S13H+S13L)/2,pressure_for_S,(S33H+S33L)/2)
plot(axial_press,S44)
xlim([0 45])
box on
hold off

figure
E3H = 1./S33H
E3L = 1./S33L

E1H = 1./S11H
E1L = 1./S11L

v13H = -1*S13H./S11H
v13L = -1*S13L./S11L
v31H = -1*S13H./S33H
v31L = -1*S13L./S33L

v12H = -1*S12H./S11H
v12L = -1*S12L./S11L


v12H(isnan(v12H))=0
v12L(isnan(v12L))=0

v13H(isnan(v13H))=0
v13L(isnan(v13L))=0


v31H(isnan(v31H))=0
v31L(isnan(v31L))=0


pgonE3 = polyshape([pressure_for_S(2:end) fliplr(pressure_for_S(2:end))],[E3H(2:end) fliplr(E3L(2:end))] )
pgonE1 = polyshape([pressure_for_S(2:end) fliplr(pressure_for_S(2:end))],[E1H(2:end) fliplr(E1L(2:end))] )

pgonv13 = polyshape([pressure_for_S(2:end) fliplr(pressure_for_S(2:end))],[v13H(2:end) fliplr(v13L(2:end))] )
pgonv12 = polyshape([pressure_for_S(2:end) fliplr(pressure_for_S(2:end))],[v12H(2:end) fliplr(v12L(2:end))] )
pgonv31 = polyshape([pressure_for_S(2:end) fliplr(pressure_for_S(2:end))],[v31H(2:end) fliplr(v31L(2:end))] )

figure
subplot(211)
hold on
plot(pgonE3)
plot(pgonE1)
xlim([5 45])
grid on
box on

subplot(212)
hold on
plot(pgonv13)
plot(pgonv12)
plot(pgonv31)
xlim([5 45])
box on
grid on
hold off






function [ep1, ep3] = V_strain_fix(ep1,ep_theta,ep3, theta)
% this function could benifit from the usage of symolic math package(e.g,
% use command: digits(50))
% following ep_theta = epV*cos(theta)^2 + epH*sin(theta)^2

if isnan(ep1)
    % with epV, ep_theta and theta to solve for epH
    ep1 = (ep_theta - ep3*(sin(theta)^2))/(cos(theta)^2);
    
elseif isnan(ep3)
    
    ep3 = (ep_theta - ep1*(cos(theta)^2))/(sin(theta)^2);
end
    
end

function [e2, e3,omega] = H_strain_gages_inversion(eH1, eH2, eH3, omega, tau1, tau2)
% this function could benifit from the usage of symolic math package(e.g,
% use command: digits(50))
% see https://www.mathworks.com/matlabcentral/answers/18991-passing-arguments-into-fsolve-without-using-globals
if isnan(omega)
% no knowledge of omega
tau = [tau1,tau2]
eH = [eH1, eH2, eH3]
f = @(x) radial_strain(x,tau,eH)
x0 = [pi/4,eH1,eH2];
options = optimoptions('fsolve','FiniteDifferenceType','central','StepTolerance',1e-20,'FunctionTolerance',1e-20,'DiffMaxChange',1e-6,'Algorithm','trust-region')
x = fsolve(f,x0,options)
omega = x(1)
e2 = x(2);
e3 = x(3);
elseif isnan(eH1)
    
elseif isnan(eH2)
    
elseif isnan(eH3)
    
end

end


function F = radial_strain(x,tau,eH)

F(1) = x(2)*(cos(x(1))^2) +x(3)*(sin(x(1))^2) - eH(1);
F(2) = x(2)*(cos(x(1)+tau(1))^2) +x(3)*(sin(x(1)+tau(1))^2) - eH(2);
F(3) = x(2)*(cos(x(1)+tau(1)+tau(2))^2) +x(3)*(sin(x(1)+tau(1)+tau(2))^2) - eH(3);

end


