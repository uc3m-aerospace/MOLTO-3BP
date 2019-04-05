%%%% State Transition Matrix for CR3BP Computation %%%%

function dqdt = DiffCorrection(t,q)

% This function will solve the following ODEs system
%                x'  = f(x)                   6x1
%               Phi' = Df(x)*Phi(t,t0)        6x6

% Control -> Choice case 

opt=2; % 1 -> Sun - Earth+Moon System
       % 2 -> Earth - Moon Syst

% State vector

% Dyn var       
x  = q(1);         
y  = q(2);         
z  = q(3);              
vx = q(4);         
vy = q(5);         
vz = q(6);         

% STM (State Transition Matrix)
Phi=zeros(6);    % Matrix initialization
Vec  = q(7:42);  % Aux vector 
for i=1:6
    for j=1:6
        k=j+6*(i-1);
        Phi(i,j)=Vec(k);
    end
end
 
%% Data

mE=5.9722e24;
mM=7.342e22;
mS=1.98845114e30;

if opt==1    % Sun-Earth
    mu=(mE+mM)/(mE+mM+mS);
else 
    mu=0.012150585609624;%1.215e-2;
end

mu1 = 1-mu;
mu2 = mu;

%% Matrix Computation

% % x' = f(x)

Ux = (mu2*(2*mu1 - 2*x))/2 - (mu1*(2*mu2 + 2*x))/2 ...
    - (mu2*(2*mu1 - 2*x))/(2*((mu1 - x)^2 + y^2 + z^2)^(3/2)) ...
    + (mu1*(2*mu2 + 2*x))/(2*((mu2 + x)^2 + y^2 + z^2)^(3/2));
Uy = (mu1*y)/((mu2 + x)^2 + y^2 + z^2)^(3/2) - mu2*y - mu1*y ...
    + (mu2*y)/((mu1 - x)^2 + y^2 + z^2)^(3/2);
Uz = (mu1*z)/((mu2 + x)^2 + y^2 + z^2)^(3/2) + ...
    (mu2*z)/((mu1 - x)^2 + y^2 + z^2)^(3/2);

% % Phi' = Df(x)*Phi(t,t0)

% Matrix U: Uxx, Uyy, Uzz .... (Sign - is still involved in this matrix)
    % U is equivalent to -U matrix from Koon ref. (P161)

m11 = (3*mu1*(2*mu2 + 2*x)^2)/(4*((mu2 + x)^2 + y^2 + z^2)^(5/2)) ...
    - mu1/((mu2 + x)^2 + y^2 + z^2)^(3/2) ...
    - mu2/((mu1 - x)^2 + y^2 + z^2)^(3/2) ...
    + (3*mu2*(2*mu1 - 2*x)^2)/(4*((mu1 - x)^2 + y^2 + z^2)^(5/2)) + 1;

m12 = (3*mu1*y*(2*mu2 + 2*x))/(2*((mu2 + x)^2 + y^2 + z^2)^(5/2)) ...
    - (3*mu2*y*(2*mu1 - 2*x))/(2*((mu1 - x)^2 + y^2 + z^2)^(5/2));

m13 = (3*mu1*z*(2*mu2 + 2*x))/(2*((mu2 + x)^2 + y^2 + z^2)^(5/2)) ...
    - (3*mu2*z*(2*mu1 - 2*x))/(2*((mu1 - x)^2 + y^2 + z^2)^(5/2));

m21 =(3*mu1*y*(2*mu2 + 2*x))/(2*((mu2 + x)^2 + y^2 + z^2)^(5/2)) ...
    - (3*mu2*y*(2*mu1 - 2*x))/(2*((mu1 - x)^2 + y^2 + z^2)^(5/2));

m22 = (3*mu2*y^2)/((mu1 - x)^2 + y^2 + z^2)^(5/2) ...
    - mu1/((mu2 + x)^2 + y^2 + z^2)^(3/2) ...
    - mu2/((mu1 - x)^2 + y^2 + z^2)^(3/2) ...
    + (3*mu1*y^2)/((mu2 + x)^2 + y^2 + z^2)^(5/2) + 1;

m23 = (3*mu2*y*z)/((mu1 - x)^2 + y^2 + z^2)^(5/2) ...
    + (3*mu1*y*z)/((mu2 + x)^2 + y^2 + z^2)^(5/2);

m31 = (3*mu1*z*(2*mu2 + 2*x))/(2*((mu2 + x)^2 + y^2 + z^2)^(5/2)) ...
    - (3*mu2*z*(2*mu1 - 2*x))/(2*((mu1 - x)^2 + y^2 + z^2)^(5/2));

m32 = (3*mu2*y*z)/((mu1 - x)^2 + y^2 + z^2)^(5/2) ...
    + (3*mu1*y*z)/((mu2 + x)^2 + y^2 + z^2)^(5/2);

m33 = (3*mu2*z^2)/((mu1 - x)^2 + y^2 + z^2)^(5/2) ...
    - mu1/((mu2 + x)^2 + y^2 + z^2)^(3/2) ...
    - mu2/((mu1 - x)^2 + y^2 + z^2)^(3/2) ...
    + (3*mu1*z^2)/((mu2 + x)^2 + y^2 + z^2)^(5/2);

U=zeros(3);   U(1,1) = m11;     U(1,2) = m12;   U(1,3) = m13;
              U(2,1) = m21;     U(2,2) = m22;   U(2,3) = m23;
              U(3,1) = m31;     U(3,2) = m32;   U(3,3) = m33;
              
% Omega Matrix + Identity matrix (I3) + zeros 3x3 matrix (Z3)
    % Equivalent to 2*Omega matrix in Koon ref.(P161)
    
Omega = [ 0   2   0;
         -2   0   0;
          0   0   0];
     
I3 = eye(3);          Z3 = zeros(3);
    
% Df(x) Matrix

Df = [Z3   I3;      % 6x6 Matrix
      U   Omega];
  
%% Ode Work

% x' = f(x)

xdot = vx;
ydot = vy;
zdot = vz;
xdotdot =  2*vy-Ux;
ydotdot = -2*vx-Uy;
zdotdot = -Uz;

% Phi' = Df(x)*Phi(t,t0)

Phidot=Df*Phi;

%% State Vector Derivative

dqdt=zeros(size(q));
dqdt(1)=xdot;
dqdt(2)=ydot;
dqdt(3)=zdot;
dqdt(4)=xdotdot;
dqdt(5)=ydotdot;
dqdt(6)=zdotdot;

PhidotAux=size(36,1); % Aux vector
for i=1:6
    for j=1:6
        k=j+6*(i-1);
        PhidotAux(k)=Phidot(i,j);
    end
end
dqdt(7:42)=PhidotAux;