%%%%% HALO ORBITS GENERATION %%%%%

%%%% 1st approx - Linear equations solution

%%% Choice case 

opt=1; % 1 -> Sun - Earth+Moon System
       % 2 -> Earth - Moon Syst
       
m=1;   % 1 -> Nothern orbit   A_z > 0
       % 3 -> Southern orbit  A_z < 0

phi=0; % rad  <-- Select initial value
Az=110e3; % km  <-- Select initial value
Ax=206e3; % km 

t=linspace(0,100,1e4); % s <-- Propagation time
       
% % DATA 

mE=5.9722e24;
mM=7.342e22;
mS=1.98845114e30;

if opt==1    % Sun-Earth
    mu=(mE+mM)/(mE+mM+mS);%3.036e-6;
else 
    mu=1.215e-2;
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
gamma2=root1(imag(root2)==0);

    % cn coefficients calculation
cn1= @(n) 1/gamma1^3*(+1.^n*mu+(-1).^n.*((1-mu)*gamma1.^(n+1))./(1-gamma1).^(n+1)); % L1
cn2= @(n) 1/gamma2^3*(-1.^n*mu+(-1).^n.*((1-mu)*gamma2.^(n+1))./(1+gamma2).^(n+1)); % L2

if opt==1 % Check with Koon results
    disp('---- CHECKs - with Koon ------')
    fprintf('\n')
    disp(['c_2 = ', num2str(cn1(2),11)])
    disp(['c_3 = ', num2str(cn1(3),11)])
    disp(['c_4 = ', num2str(cn1(4),11)])
    fprintf('\n')
end

% Eigenvalues Calculation

c2=cn1(2); % PROVISIONAL PARA COMPROBAR CON KOON
c3=cn1(3); c4=cn1(4);

% lambda=sqrt((c2+sqrt(9*c2^2-8*c2))/2);
wp=sqrt((2-c2+sqrt(9*c2^2-8*c2))/2);
wv=sqrt(c2);
lambda=wp;

if opt==1 % Check with Koon results
    disp(['w_p = ', num2str(wp,10)])
    disp(['w_v = ', num2str(wv,11)])
    fprintf('\n')
end

% Variables - Linear solution

kappa=(wp^2+1+2*c2)/(2*wp);
kappa_C=2*lambda/(lambda^2+1-c2);

% Third-Order Richardson Expansion - Intro

d1 = 3*lambda^2/kappa*(kappa*(6*lambda^2-1)-2*lambda);

a21 = (3*c3*(kappa^2-2))/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
a23 = -3*c3*lambda/(4*kappa*d1)*(3*kappa^3*lambda-6*kappa*(kappa-lambda)+4);
a24 = -3*c3*lambda/(4*kappa*d1)*(2+3*kappa*lambda);

b21 = -3*c3*lambda/(2*d1)*(3*kappa*lambda-4);
b22 = 3*c3*lambda/d1;

d21 = -c3/(2*lambda^2);

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

% Extra    - for 3rd order Richardson

d2=8*lambda^2/kappa*(kappa*(11*lambda^2-1)-2*lambda);

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

%  Extra

% Relationships

    % Amplitude Constraint
AxC=sqrt(-(incre+l2*Az^2)/l1); % km
if AxC<Ax
    disp('Amplitude Constraint is Satisfied')
else
    disp('FAILED Amplitude Constraint!!')
end

    % Phase Angle Relationship
psi=phi+m*pi/2; % rad

% EoM

x = -Ax.*cos(wp.*t+phi);
y = kappa*Ax.*sin(wp.*t+phi);
z = Az.*sin(wv.*t+psi);

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
