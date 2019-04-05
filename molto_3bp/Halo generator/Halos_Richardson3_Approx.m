%%%%% HALO ORBITS GENERATION %%%%%

%%%% 1st approx - 3rd Order Richardson Expansion

%%% Choice case 

opt=1; % 1 -> Sun - Earth+Moon System
       % 2 -> Earth - Moon Syst
       
LP=1; % 1 -> Lagrange point L1
      % 2 -> Lagrange point L2
       
m=1;   % 1 -> Nothern orbit   A_z > 0
       % 3 -> Southern orbit  A_z < 0

phi=0; % rad  <-- Select initial value
AzDim=220e3; % km  <-- Select initial value
AxDim=206e3; % km (This value will not be fixed, as Ax=f(Az))

t=linspace(0,10,1e3); % s <-- Propagation time

% % DATA 

mE=5.9722e24;
mM=7.342e22;
mS=1.98845114e30;

if opt==1    % Sun-Earth
    mu=(mE+mM)/(mE+mM+mS);%3.036e-6;
    
    % Normalization relations
                     %Intro to Lpoints
    L=1.497610041e6; %1.845287088696577e+06; %km 1.428571428571428e+06%   REVISAR!!!!!!!!!!!!!!
    Ltheory=1.496e8;
    TSE=3.155815e7; %3.2035e7; % s                         REVISAR!!!!!!!!!!!!!!
    Tconversion=TSE/(2*pi);
else 
    mu=0.012150585609624;%1.215e-2;
    L=3.84388174e5;      % REVISAR!!!!!!!!!!!!!!
    Ltheory=3.84388174e5;
    TSE=2.361e6; % s     REVISAR!!!!!!!!!!!!!!                    
    Tconversion=TSE/(2*pi);
end

% Basic Variables Calculation

rh=(mu/3)^(1/3);
gamma1ch=rh*(1-1/3*rh-1/9*rh^2-1/27*rh^3);
gamma2ch=rh*(1+1/3*rh-1/9*rh^2);

    % gamma polynomial
poly =[1 -(3-mu) 3-2*mu -mu +2*mu -mu]; %for gamma 1 --> L1
poly2=[1 +(3-mu) 3-2*mu -mu -2*mu -mu]; %for gamma 2 --> L2

root1=roots(poly);
root2=roots(poly2);
gamma1=root1(imag(root1)==0);
gamma2=root2(imag(root2)==0);

    % cn coefficients calculation
cn1= @(n) 1/gamma1^3*((+1).^n*mu+(-1).^n.*((1-mu)*gamma1.^(n+1))./(1-gamma1).^(n+1)); % L1
cn2= @(n) 1/gamma2^3*((-1).^n*mu+(-1).^n.*((1-mu)*gamma2.^(n+1))./(1+gamma2).^(n+1)); % L2

if LP==1
    c2=cn1(2); % PROVISIONAL PARA COMPROBAR CON KOON
    c3=cn1(3); c4=cn1(4);
else
    c2=cn2(2); % PROVISIONAL PARA COMPROBAR CON KOON
    c3=cn2(3); c4=cn2(4);
end

r1=Ltheory*((1-mu)-gamma1);  % REVISAR!!!!!!!!!!!!
r2=Ltheory*((1-mu)+gamma2);

if opt==1 % Check with Koon results
    disp('---- CHECKs - with Koon ------')
    fprintf('\n')
    disp(['c_2 = ', num2str(c2,11)])
    disp(['c_3 = ', num2str(c3,11)])
    disp(['c_4 = ', num2str(c4,11)])
    fprintf('\n')
end

% Eigenvalues Calculation

% lambda=sqrt((c2+sqrt(9*c2^2-8*c2))/2);
wp=sqrt((2-c2+sqrt(9*c2^2-8*c2))/2);
wv=sqrt(c2);
lambda=wp; % For halo orbits

if opt==1 % Check with Koon results
    disp(['w_p = ', num2str(wp,10)])
    disp(['w_v = ', num2str(wv,11)])
    fprintf('\n')
end

% Variables - Linear solution

kappa=(wp^2+1+2*c2)/(2*wp);
kappa_C=2*lambda/(lambda^2+1-c2);

% Third-Order Richardson Expansion 

d1 = 3*lambda^2/kappa*(kappa*(6*lambda^2-1)-2*lambda);
d2= 8*lambda^2/kappa*(kappa*(11*lambda^2-1)-2*lambda);

a21 = (3*c3*(kappa^2-2))/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
a23 = -3*c3*lambda/(4*kappa*d1)*(3*kappa^3*lambda-6*kappa*(kappa-lambda)+4);
a24 = -3*c3*lambda/(4*kappa*d1)*(2+3*kappa*lambda);

b21 = -3*c3*lambda/(2*d1)*(3*kappa*lambda-4);
b22 = 3*c3*lambda/d1;

d21 = -c3/(2*lambda^2);

a31=-9*lambda/(4*d2)*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa^2))+...
    (9*lambda^2+1-c2)/(2*d2)*(3*c3*(2*a23-kappa*b21)+c4*(2+3*kappa^2));
a32=-1/d2*(9*lambda/4*(4*c3*(kappa*a24-b22)+kappa*c4)+...
    3/2*(9*lambda^2+1-c2)*(c3*(kappa*b22+d21-2*a24)-c4));

b31=3/(8*d2)*(8*lambda*(3*c3*(kappa*b21-2*a23)-c4*(2+3*kappa^2))+...
    (9*lambda^2+1+2*c2)*(4*c3*(kappa*a23-b21)+kappa*c4*(4+kappa^2)));
b32=1/d2*(9*lambda*(c3*(kappa*b22+d21-2*a24)-c4)+...
    3/8*(9*lambda^2+1+2*c2)*(4*c3*(kappa*a24-b22)+kappa*c4));

d31=3/(64*lambda^2)*(4*c3*a24+c4);
d32=3/(64*lambda^2)*(4*c3*(a23-d21)+c4*(4+kappa^2));

s1 = (2*lambda*(lambda*(1+kappa^2)-2*kappa))^-1*...
    (3/2*c3*(2*a21*(kappa^2-2)-a23*(kappa^2+2)-2*kappa*b21)-3/8*c4*(3*kappa^4-8*kappa^2+8));
s2  = (2*lambda*(lambda*(1+kappa^2)-2*kappa))^-1*...
    (3/2*c3*(2*a22*(kappa^2-2)+a24*(kappa^2+2)+2*kappa*b22+5*d21)+3/8*c4*(12-kappa^2));


l1=-3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-kappa^2)+2*lambda^2*s1;
l2=3/2*c3*(a24-2*a22)+9/8*c4+2*lambda^2*s2;
incre=wp^2-wv^2;

if opt==1 % Check with Koon results
    disp(['s_1 = ', num2str(s1,10)])
    disp(['s_2 = ', num2str(s2,10)])
    disp(['l_1 = ', num2str(l1,10)])
    disp(['l_2 = ', num2str(l2,10)])
    disp(['Increment = ', num2str(incre,11)])
    fprintf('\n')
end

% Relationships

    % Amplitude Constraint

Ax=AxDim/L; Az=AzDim/L;
AxC=sqrt(-(incre+l2*Az^2)/l1); % km
if AxC==Ax
    disp('Amplitude Constraint is Satisfied')
else
    disp('FULFILLING Amplitude Constraint!!')   % REVISAR!!!!!!!!
    Ax=AxC;
end
    % Phase Angle Relationship
psi=phi+m*pi/2; % rad

% Added variables - Lindstedt-Poncairé Method

nu1=0;
nu2=s1*(Ax)^2+s2*(Az)^2; % REVISAR!!!!!!!!!!!
nu=1+nu1+nu2;
tau=nu.*t;

tau1=wp.*tau+phi;
deltam=2-m;

T=2*pi/(wp*nu)*Tconversion; % Period Calculation (s)
Tdays=T/3600/24;    % s to days

% EoM - 3rd Order Richardson Expansion

xad = a21*Ax^2+a22*Az^2-Ax.*cos(tau1)+...
    (a23*Ax^2-a24*Az^2).*cos(2*tau1)+(a31*Ax^3-a32*Ax*Az^2).*cos(3*tau1);
yad = kappa*Ax.*sin(tau1)+...
    (b21*Ax^2-b22*Az^2).*sin(2*tau1)+(b31*Ax^3-b32*Ax*Az^2)*sin(3*tau1);
zad = deltam*Az.*cos(tau1)+...
    deltam*d21*Ax*Az.*(cos(2*tau1)-3)+deltam*(d32*Az*Ax^2-d31*Az^3).*cos(3*tau1);

x=xad*L; y=yad*L; z=zad*L;
% % Plots 

figure()
plot3(x,y,z); xlabel('x'); ylabel('y'); zlabel('z');

figure()
subplot(1,3,1)
plot(x,y); xlabel('x'); ylabel('y');
subplot(1,3,2)
plot(x,z); xlabel('x'); ylabel('z');
subplot(1,3,3)
plot(y,z); xlabel('y'); ylabel('z');

    % 1 orbit display

tdim=t*Tconversion;
tf=min(find(tdim>=T));
figure()
subplot(1,3,1)
plot(x(1:tf),y(1:tf)); xlabel('x'); ylabel('y'); grid on; xlim([-10e5 10e5]); ylim([-10e5 10e5]); %axis equal;
hold on
plot(0,0,'bx')
subplot(1,3,2)
plot(x(1:tf),z(1:tf)); xlabel('x'); ylabel('z'); grid on; xlim([-10e5 10e5]); ylim([-10e5 10e5]); %axis equal;
hold on
plot(0,0,'bx')
subplot(1,3,3)
plot(y(1:tf),z(1:tf)); xlabel('y'); ylabel('z'); grid on; xlim([-10e5 10e5]); ylim([-10e5 10e5]); %axis equal;
hold on
plot(0,0,'bx')
if m==1
    suptitle(['Class I  Orbit  (L',num2str(LP),')'])
else
    suptitle(['Class II Orbit  (L',num2str(LP),')'])
end
% figure() % Case Intro to Lpoints
% subplot(1,3,1)
% plot(r1+x(1:tf),y(1:tf)); xlabel('x'); ylabel('y'); xlim([1.47e8 1.49e8]); grid on;%axis equal;
% hold on
% plot(r1,0,'bx')
% subplot(1,3,2)
% plot(r1+x(1:tf),z(1:tf)); xlabel('x'); ylabel('z'); xlim([1.47e8 1.49e8]); ylim([-7e5 7e5]); grid on;%axis equal;
% hold on
% plot(r1,0,'bx')
% subplot(1,3,3)
% plot(y(1:tf),z(1:tf)); xlabel('y'); ylabel('z'); xlim([-12e5 12e5]); ylim([-7e5 7e5]); grid on;
% hold on
% plot(0,0,'bx')

% Displays

fprintf('\n')
if opt==1
    disp('RCTBP: SUN-(EARTH+MOON) SYSTEM')
else
    disp('RCTBP: EARTH-MOON SYSTEM')
end
disp('-- Results ---')
disp(['Ax = ', num2str(AxC*L), ' km;   Az = ', num2str(AzDim), ' km'])
disp(['T (days)= ', num2str(Tdays,6)])
if m==1
    disp(['Case: Northern (L',num2str(LP),')'])
else
    disp(['Case: Southern (L',num2str(LP),')'])
end
fprintf('\n')