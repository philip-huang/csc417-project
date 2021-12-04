close all
clear all;
clear
%% Rigid-body bar link parameters
COM_init = [0 0]; % center of mass
w = 1; % width of rigid body
h = 5; % height of rigid body
j = 0.25; % joint hole size
d = 2; % distance of the joint holes away from COM, the joint holes on 
       % both sides are the same distance away from the COM
m = 5; % mass of the object

% Define translation frame of each rigid body
xL = [1;0];
yL = [0;1];
x = [xL yL];

%% Visualize Test
q0 = [0 0 0 0 0 2*d]'; % initial position of the system (6x1)
% (th_a, x_a, y_a, th_b, x_b, y_b) th [degrees]
% x,y position of COM of each body

% figure (2)
% hold on
% visualize(q0);
% axis equal
%% Get the Jacobian
syms p1_a p2_a p1_b p2_b real % joint coordinates in the body frame***

% gen coords for rigid body B, COM position and orientation about the COM
syms th_a p1_ac p2_ac th_b p1_bc p2_bc real 
q_a = [th_a; p1_ac; p2_ac]; 
q_b = [th_b; p1_bc; p2_bc];
q = [q_a; q_b];

syms w_a v1_ac v2_ac w_b v1_bc v2_bc real % angular velocity and velocity of the rigid body COM in the world frame 
qdot_a = [w_a; v1_ac; v2_ac]; 
qdot_b = [w_b; v1_bc; v2_bc];
qdot = [qdot_a; qdot_b];

syms m % uniform mass for all bodies

X_brac = @(x,y) [0 0 y; 0 0 -x; -y x 0]; % skew function

R = @(theta) [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];

v_a = R(th_a)*X_brac(p1_a, p2_a)'*R(th_a)'*[0;0;w_a];
v_b = R(th_b)*X_brac(p1_b, p2_b)'*R(th_b)'*[0;0;w_b];

dRdt = @(theta) [-sin(theta) cos(theta);
              -cos(theta) -sin(theta)];

v2d_a = dRdt(th_a) * [p1_a; p2_a] * w_a;
v2d_b = dRdt(th_b) * [p1_b; p2_b] * w_b;
% x_p = R(the) * p_b +  x_c
% v_p = dR/dtheta * w * p_b + v_c
%%
% get the velocity of the joint in world space
v_a = qdot_a(2:3) + v2d_a;
v_b = qdot_b(2:3) + v2d_b;

% Constraint
C = (v_a - v_b);

% Constraint Jacobian
J = [jacobian(C, qdot_a) jacobian(C, qdot_b)]; % (2x6)
%J = [jacobian(C, qdot)]; % (2x6)
Jfunc = matlabFunction(J);

%% Time Integration
Fext = [0 0 0 0 0 5]';
dt = 0.1;
time = 3;
tall = 1:dt:time;
stps = time/dt;
m = 1;
M = eye(6)*m;
p1_a_val = 0;
p2_a_val = d;
p1_b_val = 0;
p2_b_val = -d;

A = J*inv(M)*J';
Afunc = matlabFunction(A);

q_val = q0; % initial position of the system
qdot_val = [0 0 0 0 0 0]'; % initial velocity of the system

allConstrF = zeros(6,stps);
extF = zeros(6,stps);

%% Run 
vid = VideoWriter('video.avi');
open(vid);

figure(1)
hold on
xlim([-15 15])
ylim([-5 20])
grid on
i = 1;
for t=tall
    t
    Jval = Jfunc(p1_a_val,p1_b_val,p2_a_val,p2_b_val,q_val(1),q_val(4));
    visualize(q_val);
    bval = -(Jval*inv(M)*Fext);
    Aval = Afunc(p1_a_val,p1_b_val,p2_a_val,p2_b_val,q_val(1),q_val(4));
    lambda = Aval\bval; %inv(Aval)*bval;
    qddot_val = inv(M)*Jval'*lambda + inv(M)*Fext; % (6x1) update the velocity
    allConstrF(:,i) = inv(M)*Jval'*lambda;
    extF(:,i) = inv(M)*Fext;
    [q_val, qdot_val] = forwardeuler(q_val, qdot_val, qddot_val, dt);
    
    frame = getframe(gcf);
    writeVideo(vid,frame);
    cla
    i=i+1;
end
close(vid);

%%
% Plot the constraint forces 

figure(3)

subplot(3,1,1)
hold on
title("Forces in x")
plot(allConstrF(2, :), 'r')
plot(allConstrF(5, :), 'b')
plot(extF(2, :), 'm')
plot(extF(5, :), '--g')
legend('A constraint', 'B constraint', 'A external', 'B external')

subplot(3,1,2)
hold on
title("Forces in y")
plot(allConstrF(3, :), 'r')
plot(allConstrF(6, :), 'b')
plot(extF(3, :), 'm')
plot(extF(6, :), '--g')
legend('A constraint', 'B constraint', 'A external', 'B external')

subplot(3,1,3)
hold on
title("Torques")
plot(allConstrF(1, :), 'r')
plot(allConstrF(4, :), 'b')
plot(extF(1, :), 'm')
plot(extF(4, :), '--g')
legend('A constraint', 'B constraint', 'A external', 'B external')

