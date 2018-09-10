%Da condutividade até a constante de acoplamento
clc
clear all
%Contanted fisicas
planck = 6.62606957e-34;     % Constante de Planck [J*s]
kB = 1.3806488e-23;          % Constante de Boltzmann [J/K]
h_ = planck/(2*pi);          % Constante de Planck reduzida [J*s]
e0 = 8.854187817620e-12;     % permittivity in vacuum, [F/m]
u0 = 1.2566e-6;              % Permiabilidade do vácuo [H/m]
c = 299792458;               % Velocidade da luz no vacuo [m/s]
e = 1.60217e-19;             % Carga elementar do eletron [C]
No = 376.73;                 % Impedância intrínseca do espaço livre [Homn]
%Constantes Materiais
espd = 2.1;                  %Constante dieletrica do meio de silica
d = 50e-9;                   %Distancia entre as folhas de Grafeno [m]
T = 300;                     %Temperatura [K]
lam = 10e-6;                 %Comprimento de onda [m]
ko = 2*pi/lam;               %Vetor de onda no vacuo
w=2*pi*(c/lam);              % Frequência ângular [rad/s]
delta = 1e-9;                % Epesura da folha de Grafeno [m]
% La = 200                   % Distancia orizontal entre as folhas [m]
t = 0.5e-12;                 % Tempo de relaxação [ps]
Vf = 1e6;                    % Velocidade de Fermi [m/s]
uC = 0.2*e;                  % Potencial quimico [eV]

%%%%%%%%%%%%%%%%%%%%%%%%%% Condutividade [S]  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================%

A = (2*uC)-(w+1i*t^-1)*h_;
B = (2*uC)+(w+1i*t^-1)*h_;

% G = ((h_*(w+1i*t^-1))-2*uC)/((h_*(w+1i*t^-1))+2*uC);
H = 2*ko*pi*h_^2/(No*e^2*uC);                  % Usado para determinar kSPP
% cond_inter = (e^2/(4*h_))*(1+(1i/pi)*log(G));            %A. Wirth L. Jr.

% C = uC/(kB*T);
% cond_inter = (-1i*e^2/(4*pi*h_))*log(A/B);
%   
% cond_intra = (1i*e^2*kB*T/(pi*h_^2*w-pi*h_^2*2i*t^-1))*(C+2*log(exp(-C)+1));


cond_intra = (1i*e^2*uC)/((pi*h_^2)*(w+1i*t^-1));         %Kelvin J. A. Ooi

cond_inter = (1i*e^2/(4*pi*h_))*log(A/B);                 %Kelvin J. A. Ooi
  

sig = cond_intra + cond_inter;

sig3D = sig/delta;

Ep = 1+(1i*sig/(w*e0*delta))
Ep1 = 1+(1i*sig*No/(ko*delta));


%%%%%%%%%%%%%%%% Constantes de acoplamento [rad]  %%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================%

% d = 0:0.01e-9:80e-9;

% for k=1:length(d)
    
Kspp = ko*sqrt(espd-(2*espd/(No*sig))^2);
kSPP = H*w+(H/t)*1i;
Kp = (Kspp^2-espd*ko^2)^(1/2);

up = exp(-Kp*d);
D = 2i*espd*ko/(No*sig);
E = ((1-up)*Kspp)/(Kp+up*Kspp*d);
F = ((1+up)*Kspp)/(Kp-up*Kspp*d);

Beta_even = Kspp + ((D-Kp*(1-up))/E);

Beta_odd = Kspp + ((D-Kp*(1+up))/F);

Cg = abs((Beta_odd - Beta_even)/2)
% end

Lc1 = pi/(2*sqrt(2)*abs(Cg))
Lc2 = pi/(2*abs(Cg));

%%%%%%%%%%%%%%%%%%% Transmitância e refletância  %%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================%
