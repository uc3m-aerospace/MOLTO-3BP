%%% Main 2 - Halo Orbit Numerical Computation - x0 is kept constant

% Numerical Computation characteristics

nmax=50; % Max. number of iteration
tol=1e-15; % Tolerance error

% Orbit guess (Establish max for ode work)
tf=3;

%% IC guess

% Use Results from Main1_Halos_R3O.m

%--- Initial Guess Display ---
x0  = 1.142198291366583	;%0.98887628040431608;
y0  = 0;
z0  = -0.1599;%0.0009208314099650594;
vx0 = 3.9358748443008021e-05;
vy0 = -0.223;%0.0089026150891582875;
vz0 = -1.5451719395728401e-05;

% Recall!!! y0,vx0,vz0 = 0 !!!!!
y0=0;    vx0=0;      vz0=0;

%% Differential Correction

for i=1:nmax
    
    % IC vector preparation
    
    % x' = f(x)    IC
    q0=zeros(1,42); %Initial Conditions
    q0(1)=x0;           q0(4)=vx0;
    q0(2)=y0;           q0(5)=vy0;
    q0(3)=z0;           q0(6)=vz0;
    
    % Phi' = Df(x)*Phi(t,t0)     IC
    Phi0=eye(6);
    Phi0Vec = MattoVec(Phi0); % Aux vector
    q0(7:42)=Phi0Vec;
    
    
    % Ode45 - State Transition Matrix Computation
    
    % Event to stop ode45, x-z plane cross (at T/2)
    Opt = odeset('Events', @EventHalo);
    
    [t,q]= ode45(@DiffCorrection,[0 tf],q0,Opt);
    
    % Extracting solution
    
    xfvec   = q(end,1:6);
    xbfvec  = q(end-1,1:6);
    
    Phifvec = q(end,7:42);
    Phifmat = VectoMat6(Phifvec); % to matrix form
    
    % Desired values at tf: vxf,vzf = 0   (yf=0 already with stop event)
    % Desired values             % Final values
    vxfdes=0;  vzfdes=0;       vxf=xfvec(4);  vzf=xfvec(6);
    
    ydot=xfvec(5); %NEW
    xdotdot= (vxf-xbfvec(4))/(t(end)-t(end-1)); %NEW
    zdotdot= (vzf-xbfvec(6))/(t(end)-t(end-1)); %NEW
    
    % Delta x
    dvx=vxfdes-vxf;     dvz=vzfdes-vzf;
    B=[dvx;dvz];
    D=[xdotdot;zdotdot]; %NEW
    E=[Phifmat(2,3) Phifmat(2,5)]; %NEW
    
    % Check of IC
    
    err1 = abs(dvx);      err2 = abs(dvz);
    
    if (err1<=tol) && (err2<=tol)
        break
    else
        
        % Update IC --- Ax=B
        
        A=zeros(2,2);
        A(1,1) = Phifmat(4,3);   A(1,2) = Phifmat(4,5);
        A(2,1) = Phifmat(6,3);   A(2,2) = Phifmat(6,5);
        C=A-1/ydot*D*E; %NEW
        dxvec0=C\B; % Solve inverting A
        dz0=dxvec0(1);
        dvy0=dxvec0(2);
        
        z0  =  z0 + dz0;
        vy0 = vy0 + dvy0;
    end
    
end
%% Solution
 
disp('--- Halo Generator: Numerical Computation ---')
fprintf('\n')
if (err1<=tol) && (err2<=tol)
    disp(['Nº of iteration to converge: ', num2str(i)])
    fprintf('\n')
    disp('--- Solution ---')
    disp(['x0  = ', num2str(x0,20),';'])
    disp(['y0  = ', num2str(y0,20),';'])
    disp(['z0  = ', num2str(z0,20),';'])
    disp(['vx0 = ', num2str(vx0,20),';'])
    disp(['vy0 = ', num2str(vy0,20),';'])
    disp(['vz0 = ', num2str(vz0,20),';'])
    fprintf('\n')
    disp('--- Orbit Period ---')
    disp(['T/2 = ', num2str(t(end),20)])
    disp([' T  = ', num2str(t(end)*2,20)])
else
    disp('The program has not converged!')
    disp(['err1  = ', num2str(err1)])
    disp(['err2  = ', num2str(err2)])
    fprintf('\n')
    disp(['Nº of iterations done: ', num2str(i)])
    disp(['Tolerance: ', num2str(tol)])
    disp('Try modifying the initial guess IC ...')
    disp('  ...or modifying the number of iterations and/or the tolerance')
end