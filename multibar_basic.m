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
qtest = [0 0 0 0 0 2*d]';
%q0 = [-pi/2 0 0 0 0 2*d]'; % initial position of the system (6x1)
% (th_a, x_a, y_a, th_b, x_b, y_b) th [degrees]
% x,y position of COM of each body

% figure (2)
% hold on
% visualize(q0);
% axis equal

figure(3)
hold on
xlim([-15 15])
ylim([-5 20])
grid on
[allCOM, allBars, allax] = visualizeBar(qtest,w,h,j,d,m);
plot(allBars);
plot(allCOM(1,:),allCOM(2,:),'b o')
for i = 1:length(allCOM)
    COM = allCOM(:,i);
    ax = allax(:,:,i);
    quiver(COM(1), COM(2), ax(1,1), ax(1,2),'color',[1 0 0])
    quiver(COM(1), COM(2), ax(2,1), ax(2,2),'color',[0 1 0])
end

labels = {'A','B'};
text(allCOM(1,:),allCOM(2,:),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
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

syms alpha_a a1_ac a2_ac alpha_b a1_bc a2_bc real % angular acc and acc of the rigid body COM in the world frame 
qddot_a = [alpha_a; a1_ac; a2_ac]; 
qddot_b = [alpha_b; a1_bc; a2_bc];
qddot = [qddot_a; qddot_b];

syms m % uniform mass for all bodies
syms ks kd %spring and damping constants

% derive Jacobian from position constriant - get the same J
R2 = @(theta) [cos(theta) -sin(theta);
    sin(theta) cos(theta)];
x_a = R2(th_a) * [p1_a; p2_a] + [p1_ac; p2_ac];
x_b = R2(th_b) * [p1_b; p2_b] + [p1_bc; p2_bc];

% C(q) = 0
% Cdot(q) = dC/dq * qdot = 0
% Cddot(q) = dC/dq * qddot + dCdot/dq * qdot = 0
Csym = x_a - x_b;
Jsym = jacobian(Csym, q);

Cdot_sym = Jsym * qdot;

% c = dCdot/dq * qdot = dJ/dq * qdot * qdot
% only depends on th_a and th_b

c_sym = jacobian(Cdot_sym, q) * qdot + ks * Csym + kd * Cdot_sym;
Cddot_sym = Jsym * qddot + jacobian(Cdot_sym, q) * qdot;

cfunc = matlabFunction(c_sym);
Cfunc = matlabFunction(Csym);
Jfunc = matlabFunction(Jsym);
Cddotfunc = matlabFunction(Cddot_sym);
Cdotfunc = matlabFunction(Cdot_sym);
%% Time Integration
Fext = [-10 0 10 10 20 0]';
dt = 0.1;
time = 5;
tall = 1:dt:time;
stps = time/dt;

m = 1;
I = m * (w * w + h * h) / 12.0;
M = eye(6)*m;
M(1, 1) = I;
M(4, 4) = I;
p1_a = 0;
p2_a = d;
p1_b = 0;
p2_b = -d;

ks = 100;
kd = 10;

Asym = Jsym*inv(M)*Jsym';
Afunc = matlabFunction(Asym);

%% Run 
q0 = [0 0 0 0 0 2*d]';
q = q0; % initial position of the system
qdot = [0 0 0 0 0 0]'; % initial velocity of the system
%           w vx vy ...
allConstrF = zeros(6,stps);
% external force only applies for an initial period of time
extF = zeros(6,stps);
for t=1:1
    extF(:, t) = Fext;
end


vid = VideoWriter('video.avi');
open(vid);

figure(1)
hold on
xlim([-15 15])
ylim([-5 20])
grid on
for i=1:size(tall, 2)
    t = tall(i);
    J = Jfunc(p1_a,p1_b,p2_a,p2_b,q(1),q(4));
    %visualize(q);
    %%%%
    [allCOM, allBars, allax] = visualizeBar(q,w,h,j,d,m);
    plot(allBars);
    plot(allCOM(1,:),allCOM(2,:),'b o')
    for k = 1:length(allCOM)
        COM = allCOM(:,k);
        ax = allax(:,:,k);
        quiver(COM(1), COM(2), ax(1,1), ax(1,2),'color',[1 0 0])
        quiver(COM(1), COM(2), ax(2,1), ax(2,2),'color',[0 1 0])
    end
    labels = {'A','B'};
    text(allCOM(1,:),allCOM(2,:),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    %%%%
    c = cfunc(kd, ks, p1_a, p1_b, p2_a, p2_b, q(2), q(5), q(3), q(6), q(1), q(4), ...
        qdot(2), qdot(5), qdot(3), qdot(6), qdot(1), qdot(4));
    b = -(J*inv(M)*extF(:, i) + c);
    A = Afunc(p1_a,p1_b,p2_a,p2_b,q(1),q(4));
    lambda = A\b; %inv(A)*b;
    qddot = inv(M)*J'*lambda + inv(M)*extF(:, i); % (6x1) update the velocity
    allConstrF(:,i) = J'*lambda;

    % debug information
    Cddot = Cddotfunc(qddot(2), qddot(5), qddot(3), qddot(6), qddot(1), qddot(4), p1_a, p1_b, p2_a, p2_b, q(1), q(4), qdot(1), qdot(4));

    [q, qdot] = forwardeuler(q, qdot, qddot, dt);
    
    % debug information
    Cdot = Cdotfunc(p1_a, p1_b, p2_a, p2_b, q(1), q(4), qdot(2), qdot(5), qdot(3), qdot(6), qdot(1), qdot(4));
    C = Cfunc(p1_a, p1_b, p2_a, p2_b, q(2), q(5), q(3), q(6), q(1), q(4));
    q;

    frame = getframe(gcf);
    writeVideo(vid,frame);
    cla
end
close(vid);

%%
% Plot the constraint forces 
visualize(q);
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
