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

%% Visualize Test (opened 4-bar linkage, to extend to with aux constrs)
% (th_1, x_1, y_1, th_0, x_0, y_0, ..., th_n, x_n, y_n) th [degrees]
% n = number of elements in the system
% x,y position of COM of each body

root = [pi/2,0,0]';
relth =  [-pi/2 -pi/2 0]; % initial position of the system (6x1)
q0 = generateInitCoords(root, relth, d);
figure(3)
hold on
xlim([-15 15])
ylim([-5 20])
grid on
[allCOM, allBars, allax] = getAllBars(q0,w,h,j,d,m);
visualizeAllBars(allCOM, allBars, allax);
%% Get the Jacobian
% system definition!
root = [0,0,0]';
relth = [0]; % initial position of the system (6x1)
q0 = generateInitCoords(root, relth, d); % initial position defined the number of bodies in the system
numBod = length(q0)/3;

% syms p1_a p2_a p1_b p2_b real 
% each bar has 2 attachment points, there for (numBod)*2 joint points
% pa = p1_a p2_a
pa = sym('pa', [numBod*2 2], 'real'); % joint coordinates in the body frame***
% pb = p1_b p2_b
pb = sym('pb', [numBod*2 2], 'real');

% syms th_a p1_ac p2_ac th_b p1_bc p2_bc real 
% q_a = [th_a; p1_ac; p2_ac]; 
% q_b = [th_b; p1_bc; p2_bc];
% q = [q_a; q_b];
q = sym('q', [numBod 3], 'real'); % gen coords for rigid body, COM position and orientation about the COM

% syms w_a v1_ac v2_ac w_b v1_bc v2_bc real 
% qdot_a = [w_a; v1_ac; v2_ac]; 
% qdot_b = [w_b; v1_bc; v2_bc];
% qdot = [qdot_a; qdot_b];
qdot = sym('qdot', [numBod 3], 'real'); % angular velocity and velocity of the rigid body COM in the world frame 

% syms alpha_a a1_ac a2_ac alpha_b a1_bc a2_bc real 
% qddot_a = [alpha_a; a1_ac; a2_ac]; 
% qddot_b = [alpha_b; a1_bc; a2_bc];
% qddot = [qddot_a; qddot_b];
qddot = sym('qddot', [numBod 3], 'real'); % angular acc and acc of the rigid body COM in the world frame 

% derive Jacobian from position constraint - get the same J
R2 = @(theta) [cos(theta) -sin(theta);
    sin(theta) cos(theta)];
% x_a = R2(th_a) * [p1_a; p2_a] + [p1_ac; p2_ac];
% x_b = R2(th_b) * [p1_b; p2_b] + [p1_bc; p2_bc];

% transform the joint coordinates in the body frame to the world frame
xa = R2(q(1))*pa(2, :)' + q(1,2:3)';
xb = R2(q(4))*pb(1, :)' + q(2,2:3)';

% C(q) = 0
% Cdot(q) = dC/dq * qdot = 0
% Cddot(q) = dC/dq * qddot + dCdot/dq * qdot = 0

%Csym = x_a - x_b;
Csym = xa - xb;

% flatten everything to vectors
pa = (reshape(pa,numBod*2*2,1));
pb = (reshape(pb,numBod*2*2,1));
q = (reshape(q,numBod*3,1));
qdot = (reshape(qdot,numBod*3,1));
qddot = (reshape(qddot,numBod*3,1));
Jsym = jacobian(Csym, q);

Cdot_sym = Jsym * qdot;

% c = dCdot/dq * qdot = dJ/dq * qdot * qdot
% only depends on th_a and th_b

m = 1; % uniform mass
ks = 100;%spring and damping constants
kd = 10;
c_sym = jacobian(Cdot_sym, q) * qdot + ks * Csym + kd * Cdot_sym;
%Cddot_sym = Jsym * qddot + jacobian(Cdot_sym, q) * qdot;

%Cfunc = matlabFunction(Csym);
%cfunc = matlabFunction(c_sym);
%Jfunc = matlabFunction(Jsym);

cfuncmod = matlabFunction(c_sym, 'Vars', {pa, pb, q, qdot});
Jfuncmod = matlabFunction(Jsym, 'Vars', {pa, pb, q});

% Cddotfunc = matlabFunction(Cddot_sym);
% Cdotfunc = matlabFunction(Cdot_sym);

%% Time Integration


Fext = [-10 0 10 10 20 0]';
dt = 0.1;
time = 5;
tall = 1:dt:time;
stps = time/dt;

I = m * (w * w + h * h) / 12.0;
M = eye(6)*m;
M(1, 1) = I;
M(4, 4) = I;

Asym = Jsym*inv(M)*Jsym';
Afunc = matlabFunction(Asym, 'Vars', {pa, pb, q});
% p1_a = 0;
% p2_a = d;
% p1_b = 0;
% p2_b = -d;
pa = [0 -d 0 d 0 -d 0 d]';
pb = [0 -d 0 d 0 -d 0 d]';

%% Run 
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
    %J = Jfunc(p1_a,p1_b,p2_a,p2_b,q(1),q(4));
    J = Jfuncmod(pa, pb, q);

    [allCOM, allBars, allax] = getAllBars(q,w,h,j,d,m);
    visualizeAllBars(allCOM, allBars, allax);

    %c = cfunc(p1_a, p1_b, p2_a, p2_b, q(2), q(5), q(3), q(6), q(1), q(4), ...
    %    qdot(2), qdot(5), qdot(3), qdot(6), qdot(1), qdot(4));
    c = cfuncmod(pa, pb, q, qdot);

    b = -(J*inv(M)*extF(:, i) + c);
    A = Afunc(pa, pb, q);
    lambda = A\b; %inv(A)*b;
    qddot = inv(M)*J'*lambda + inv(M)*extF(:, i); % (6x1) update the velocity
    allConstrF(:,i) = J'*lambda;

    % debug information
    %Cddot = Cddotfunc(qddot(2), qddot(5), qddot(3), qddot(6), qddot(1), qddot(4), p1_a, p1_b, p2_a, p2_b, q(1), q(4), qdot(1), qdot(4));

    [q, qdot] = forwardeuler(q, qdot, qddot, dt);
    
    % debug information
    %Cdot = Cdotfunc(p1_a, p1_b, p2_a, p2_b, q(1), q(4), qdot(2), qdot(5), qdot(3), qdot(6), qdot(1), qdot(4));
    %C = Cfunc(p1_a, p1_b, p2_a, p2_b, q(2), q(5), q(3), q(6), q(1), q(4));
    %q;

    frame = getframe(gcf);
    writeVideo(vid,frame);
    cla
end
close(vid);