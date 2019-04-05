%%% Main 3 - Halo Orbit Numerical Computation - Orbit Propagation 

% Ode 45 - General equation around Lpoint

q0=zeros(1,6); %Initial Conditions

% Use Results from Main2_Halo_num.m

x0  = 1.1107439873579033;
y0  = 0;
z0  = 0.035680331960522345;
vx0 = 0;
vy0 = 0.20363565695006591;
vz0 = 0;

q0(1)=x0;           q0(4)=vx0;
q0(2)=y0;           q0(5)=vy0;
q0(3)=z0;           q0(6)=vz0;

tf=3.39; % Period Time INSERT
tspan=linspace(0,tf,1e5);

[t,q]=ode45(@TBP,tspan,q0);
x=q(:,1); y=q(:,2); z=q(:,3);

figure()
plot3(x,y,z); xlabel('x'); ylabel('y'); zlabel('z');

figure()
subplot(1,3,1)
plot(x,y); xlabel('x'); ylabel('y');
subplot(1,3,2)
plot(x,z); xlabel('x'); ylabel('z');
subplot(1,3,3)
plot(y,z); xlabel('y'); ylabel('z');