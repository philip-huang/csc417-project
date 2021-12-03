close all
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
R = @(theta) [cos(theta) -sin(theta);
    sin(theta) cos(theta)];

%% Visualize Test
q0 = [0 0 0 0 0 2*d]'; % initial position of the system (6x1)
% (th_a, x_a, y_a, th_b, x_b, y_b) th [degrees]
% x,y position of COM of each body

% figure (2)
% hold on
% visualize(q0);
% axis equal
%% Get the Jacobian
R = @(theta) [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];

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

v_a = R(th_a)*X_brac(p1_a, p2_a)'*R(th_a)'*[0;0;w_a];
v_b = R(th_b)*X_brac(p1_b, p2_b)'*R(th_b)'*[0;0;w_b];

% get the velocity of the joint in world space
v_a = qdot_a(2:3, :);% + v_a(1:2, :);
v_b = qdot_b(2:3, :);% + v_b(1:2, :);

% Constraint
C = v_a - v_b;

% Constraint Jacobian
%J = [jacobian(C, qdot_a) jacobian(C, qdot_b)]; % (2x6)
J = [jacobian(C, qdot)]; % (2x6)
Jfunc = matlabFunction(J);

%% Time Integration
dt = 0.5;
tall = 1:dt:3*2*pi;
m = 1;
M = eye(6)*m;
p1_a = 0;
p2_a = d;
p1_b = 0;
p2_b = -d;

A = J*inv(M)*J';
Afunc = matlabFunction(A);
dt = 0.1;

q = q0; % initial position of the system
qdot = [0 0 0 0 0 0]'; % initial velocity of the system

vid = VideoWriter('video.avi');
open(vid);

figure(1)
hold on
xlim([-15 15])
ylim([-5 20])
grid on
for t=tall
    Fext = [0 0 0 0 -5 0]';
    t
    Jval = Jfunc();%(p1_a,p1_b,p2_a,p2_b,q(1),q(3));
    visualize(q);
    bval = -(Jval*inv(M)*Fext);
    Aval = Afunc();%p1_a,p1_b,p2_a,p2_b,q(1),q(3));
    lambda = Aval\bval;
    qddot = inv(M)*Jval'*lambda + inv(M)*Fext % (6x1) update the velocity
    qdot
    [q, qdot] = forwardeuler(q, qdot, qddot, dt);
    
    frame = getframe(gcf);
    writeVideo(vid,frame);
    cla
end
close(vid);
