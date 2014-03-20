% Solving Navier Stokes Equation for Lid driven cavity problem
clc;
clear;

%  Material Properties; Non-Dimensionalized
L = 1; 
H = 1; 
imax = 32;
jmax = 32;
re = input('Enter the Reynolds Number = '); % Reynolds Number
rho = 1; 
mu = 1/re; % Dynamic Viscosity

u_top = 1;
dx = L/(imax-2);
dy = H/(jmax-2);
eps = 1e-6; % Convergence Criteria
check = 1;
scale = input('Enter the dt scale down ratio = ');

tic

% Initialization of parameters
p = zeros(jmax,imax)+1; % Pressure Grid
u = zeros(jmax,imax-1); % U velocity
v = zeros(jmax-1,imax); % V Velocity
uy_advect = zeros(jmax-1,imax-3);
vx_advect = zeros(jmax-3,imax-1);

% Setting given initial conditions specific to problem
u(1,:) = u_top;

% Time step determination specific to problem
dt = dx/u_top; % CFL Criteria
dt = dt/scale; % To be on the safer side

u_old = u; v_old = v;

% Constants used
constant_a = dt/(rho*dx*dy);
aE = (rho*dt*dy)/dx; aW = aE;
aN = (rho*dt*dx)/dy; aS = aN;
aP = aE+aW+aN+aS;
rr = 1;
rr(re==400) = 2;


while check > eps
    %% SOLVING X MOMENTUM EQUATION
    
    % MASS FLOW RATE
    ux_mass = rho * (u(2:jmax-1,1:imax-2)+u(2:jmax-1,2:imax-1)) * 0.5;
    uy_mass = rho * (v(:,2:imax-2)+v(:,3:imax-1)) * 0.5;
    
    % ADVECTION CALCULATION
    % X direction advection
    ux_advect = (max(ux_mass,0).*u(2:jmax-1,1:imax-2)) + ...
        (min(ux_mass,0).*u(2:jmax-1,2:imax-1));
    % Y direction advection
    uy_advect(2:jmax-2,:) = (max(uy_mass(2:jmax-2,:),0).*u(3:jmax-1,2:imax-2)) +...
        (min(uy_mass(2:jmax-2,:),0).*u(2:jmax-2,2:imax-2));
    % Y direction boundary advection
    uy_advect(1,:) = uy_mass(1,:).*u(1,2:imax-2);
    uy_advect(jmax-1,:) = uy_mass(jmax-1,:).*u(jmax,2:imax-2);
    % Net advection calculation
    u_advect = (ux_advect(:,2:imax-2)-ux_advect(:,1:imax-3))*dy + ...
        (uy_advect(1:jmax-2,:)-uy_advect(2:jmax-1,:))*dx;
    
    % DIFFUSION CALCULATION
    % X direction diffusion
    ux_diffuse = mu * (u(2:jmax-1,2:imax-1)-u(2:jmax-1,1:imax-2)) * (1/dx);
    % Y direction diffusion
    uy_diffuse = mu * (u(1:jmax-1,2:imax-2)-u(2:jmax,2:imax-2)) * (1/dy);
    % Y direction boundary diffusion
    uy_diffuse(1,:) = uy_diffuse(1,:) * 2;
    uy_diffuse(jmax-1,:) = uy_diffuse(jmax-1,:) * 2;
    % Net diffusion calculation
    u_diffuse = (ux_diffuse(:,2:imax-2)-ux_diffuse(:,1:imax-3))*dy + ...
        (uy_diffuse(1:jmax-2,:)-uy_diffuse(2:jmax-1,:))*dx;
    
    % SOURCE CALCULATION
    u_source =    (p(2:jmax-1,2:imax-2) - p(2:jmax-1,3:imax-1)) * dy;
    
    % PREDICTION OF U VELOCITY
    u(2:jmax-1,2:imax-2) = u(2:jmax-1,2:imax-2) +...
        (constant_a*(u_diffuse-u_advect+u_source));
    
    
    %%  SOLVING Y MOMENTUM EQUATION
    
    % MASS FLOW RATE
    vx_mass = rho * (u(2:jmax-2,:)+u(3:jmax-1,:)) * 0.5;
    vy_mass = rho * (v(1:jmax-2,2:imax-1)+v(2:jmax-1,2:imax-1)) * 0.5;
    
    % ADVECTION CALCULATION
    % X direction advection
    vx_advect(:,2:imax-2) = (max(vx_mass(:,2:imax-2),0).*v(2:jmax-2,2:imax-2)) +...
        (min(vx_mass(:,2:imax-2),0).*v(2:jmax-2,3:imax-1));
    % X direction boundary advection
    vx_advect(:,1) = vx_mass(:,1).*v(2:jmax-2,1);
    vx_advect(:,imax-1) = vx_mass(:,imax-1).*v(2:jmax-2,imax);
    % Y direction advection
    vy_advect = (max(vy_mass,0).*v(2:jmax-1,2:imax-1)) + ...
        (min(vy_mass,0).*v(1:jmax-2,2:imax-1));
    % Net advection calculation
    v_advect = (vx_advect(:,2:imax-1)-vx_advect(:,1:imax-2))*dy + ...
        (vy_advect(1:jmax-3,:)-vy_advect(2:jmax-2,:))*dx;
    
    % DIFFUSION CALCULATION
    % X direction diffusion
    vx_diffuse = mu * (v(2:jmax-2,2:imax)-v(2:jmax-2,1:imax-1)) * (1/dx);
    % X direction boundary diffusion
    vx_diffuse(:,1) = vx_diffuse(:,1) * 2;
    vx_diffuse(:,imax-1) = vx_diffuse(:,imax-1) * 2;
    % Y direction diffusion
    vy_diffuse = mu * (v(1:jmax-2,2:imax-1)-v(2:jmax-1,2:imax-1)) * (1/dy);
    % Net diffusion calculation
    v_diffuse = (vx_diffuse(:,2:imax-1)-vx_diffuse(:,1:imax-2))*dy + ...
        (vy_diffuse(1:jmax-3,:)-vy_diffuse(2:jmax-2,:))*dx;
    
    % SOURCE CALCULATION
    v_source =    (p(3:jmax-1,2:imax-1) - p(2:jmax-2,2:imax-1)) * dx;
    
    % PREDICTION OF V VELOCITY
    v(2:jmax-2,2:imax-1) = v(2:jmax-2,2:imax-1) + ...
        (constant_a*(v_diffuse-v_advect+v_source));
    
    
    %% Solving Continuity Equation
    
    divergence = ((u(2:jmax-1,2:imax-1)-u(2:jmax-1,1:imax-2))*rho*dy) +...
        ((v(1:jmax-2,2:imax-1)-v(2:jmax-1,2:imax-1))*rho*dx);
    p_corrector = zeros(jmax,imax);
    
    % Gauss Siedel Iteration
    while max(max(abs(divergence)))>eps
        
        % Pressure Correction
        for j = 2:jmax-1
            for i = 2:imax-1
                p_corrector(j,i) = ((aE*p_corrector(j,i+1))+(aW*p_corrector(j,i-1))+...
                    (aN*p_corrector(j-1,i))+(aS*p_corrector(j+1,i))-divergence(j-1,i-1)) * (1/aP);
            end
        end
        
        % Pressure boundary correction
        p_corrector(2:jmax-1,1) = p_corrector(2:jmax-1,2); 
        p_corrector(2:jmax-1,imax) = p_corrector(2:jmax-1,imax-1);
        p_corrector(1,2:imax-1) = p_corrector(2,2:imax-1); 
        p_corrector(jmax,2:imax-1) = p_corrector(jmax-1,2:imax-1);
        p = p + p_corrector;
        
        % Correction of U predicted & V predicted
        u(2:jmax-1,2:imax-2) = u(2:jmax-1,2:imax-2) + ...
            (constant_a*dy*(p_corrector(2:jmax-1,2:imax-2)-p_corrector(2:jmax-1,3:imax-1)));
        v(2:jmax-2,2:imax-1) = v(2:jmax-2,2:imax-1) + ...
            (constant_a*dx*(p_corrector(3:jmax-1,2:imax-1)-p_corrector(2:jmax-2,2:imax-1)));
        
        % Divergence Calculation
        divergence = ((u(2:jmax-1,2:imax-1)-u(2:jmax-1,1:imax-2))*rho*dy) + ...
            ((v(1:jmax-2,2:imax-1)-v(2:jmax-1,2:imax-1))*rho*dx);
    end
    
    % Convergence Criteria
    check = sqrt(mean(mean((u-u_old).^2)))+sqrt(mean(mean((v-v_old).^2)));
    
    disp(check);
    u_old = u; v_old = v;
end
toc

%% Tracing velocity from staggered to normal grid

U = zeros(jmax,imax);
V = zeros(jmax,imax);
% Velocity approximation at corners
U(1:jmax-1:jmax,1:imax-1:imax) = ((2/3)*u(2:jmax-3:jmax-1,1:imax-2:imax-1)) + ((1/3)*u(1:jmax-1:jmax,2:imax-4:imax-2)); 
V(1:jmax-1:jmax,1:imax-1:imax) = (1/3)*v(2:jmax-4:jmax-2,1:imax-1:imax) + (2/3)*v(1:jmax-2:jmax-1,2:imax-3:imax-1);
% Interior points
U(2:jmax-1,2:imax-1) = (u(2:jmax-1,1:imax-2)+u(2:jmax-1,2:imax-1))*0.5;
V(2:jmax-1,2:imax-1) = (v(1:jmax-2,2:imax-1)+v(2:jmax-1,2:imax-1))*0.5;
% Boundary points
U(1,2:imax-1) = u_top;

x = [0,linspace(dx/2,L-dx/2,imax-2),L];
y = [H,linspace(H-dy/2,dy/2,jmax-2),0];

% Velocity Vector plot
quiver(x(2:imax-1),y(2:jmax-1),U(2:jmax-1,2:imax-1),V(2:jmax-1,2:imax-1),2);

center_u = (U(:,imax/2)+U(:,(imax/2)+1))*0.5;
center_v = (V(jmax/2,:)+V((jmax/2)+1,:))*0.5;
% Ghia Results
ypA = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0];
ucA(1,:) = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.2109 -0.15662 -0.1015 -0.06434 -0.04775 -0.04192 -0.03717 0];
ucA(2,:) = [1 0.75837 0.68439 0.61756 0.55892 0.29093 0.16256 0.02135 -0.11477 -0.17119 -0.32726 -0.24299 -0.14612 -0.10338 -0.09266 -0.08186 0];
vcA(1,:) = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.100091 0.09233 0];
vcA(2,:) = [0 -0.12146 -0.15663 -0.19254 -0.22847 -0.23827 -0.44993 -0.38598 0.05186 0.30174 0.30203 0.28124 0.22965 0.2092 0.19713 0.1836 0];
xpA = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];

figure;
hold all
plot(center_u,y);
plot(ucA(rr,:),ypA,'*-');
legend('Present','Published');
title('U-Velocity along vertical centerline');
figure;
hold all
plot(x,center_v);
plot(xpA,vcA(rr,:),'*-');
legend('Present','Published')
title('V-Velocity along horizontal centerline');

figure;
contourf(x,y,U,20);
figure;
contourf(x,y,V,20);