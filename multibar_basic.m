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
m = 1; % mass of the objects

ks = 100;%spring and damping constants
kd = 10;

% Define translation frame of each rigid body
xL = [1;0];
yL = [0;1];
x = [xL yL];

%% Visualize Test (opened 4-bar linkage, to extend to with aux constrs)
% (th_1, x_1, y_1, th_0, x_0, y_0, ..., th_n, x_n, y_n) th [degrees]
% n = number of elements in the system
% x,y position of COM of each body

% root = [pi/2,0,0]';
% relth =  [-pi/2 -pi/2 0]; % initial position of the system (6x1)
% q0 = generateInitCoords(root, relth, d);
% figure(3)
% hold on
% xlim([-15 15])
% ylim([-5 20])
% grid on
% [allCOM, allBars, allax] = getAllBars(q0,w,h,j,d,m);
% visualizeAllBars(allCOM, allBars, allax);

%% Get the Jacobian
% system definition!

% % 2-bar case
% root = [0,0,0]';
% relth = [0]; % initial position of the system (6x1)
% q0 = generateInitCoords(root, relth, d); % initial position defined the number of bodies in the system
% numBod = length(q0)/3;

% 4-bar case
root = [pi/2,0,0]';
relth =  [-pi/2 -pi/2 0 0]; % initial position of the system (6x1)
q0 = generateInitCoords(root, relth, d);
numBod = length(q0)/3;

% joint coordinates in the body frame, one constraint between each joint
% NOTE: Add more joint coordinates as you add more constraints!
pa = sym('pa', [(numBod-1) 2], 'real');
pb = sym('pb', [(numBod-1) 2], 'real');

q = sym('q', [numBod 3], 'real'); % gen coords for rigid body, COM position and orientation about the COM
qdot = sym('qdot', [numBod 3], 'real'); % angular velocity and velocity of the rigid body COM in the world frame 
qddot = sym('qddot', [numBod 3], 'real'); % angular acc and acc of the rigid body COM in the world frame 

R2 = @(theta) [cos(theta) -sin(theta);
    sin(theta) cos(theta)];

% transform the joint coordinates in the body frame to the world frame
xa = [];
xb = [];
for c = 1:(numBod-1) % for every constraint, body i-1 and body i
    xa = cat(1,xa,R2(q(c,1))*pa(1,:)' + q(c,2:3)');
    xb = cat(1,xb,R2(q(c+1,1))*pb(1,:)' + q(c+1,2:3)');
end

% C(q) = 0
% Cdot(q) = dC/dq * qdot = 0
% Cddot(q) = dC/dq * qddot + dCdot/dq * qdot = 0
Csym = xa - xb;

% flatten everything to vectors
pa = reshape(pa',(numBod-1)*2,1);
pb = reshape(pb',(numBod-1)*2,1);
q = reshape(q',numBod*3,1);
qdot = reshape(qdot',numBod*3,1);
qddot = reshape(qddot',numBod*3,1);
Jsym = jacobian(Csym', q);

Cdot_sym = Jsym * qdot;

c_sym = jacobian(Cdot_sym, q) * qdot + ks * Csym + kd * Cdot_sym;
cfuncmod = matlabFunction(c_sym, 'Vars', {pa, pb, q, qdot});
Jfuncmod = matlabFunction(Jsym, 'Vars', {pa, pb, q});


%% Time Integration
%Fext = zeros(numBod*3,1); % apply zero external force

% Apply a tip force
Fext = zeros(numBod*3,1);
Fext(end-2:end) = [0 15 35]';

% Fext = 10*ones(numBod*3,1); % apply forces on all the joints

dt = 0.1;
time = 5;
tall = 1:dt:time;
stps = time/dt;

I = m * (w * w + h * h) / 12.0;
M = eye(numBod*3)*m;

for i = 0:numBod-1
    M(i*3+1,i*3+1) = I;
end

Asym = Jsym*inv(M)*Jsym';
Afunc = matlabFunction(Asym, 'Vars', {pa, pb, q});

% Give the positions of the joint that are being constrained
pa = zeros((numBod-1)*2*2,1);
pb = zeros((numBod-1)*2*2,1);
% NOTE: this is for a connection between each individual joint that forms a
% serial chain
for i = 0:(numBod-1)*2-1
    if(rem(i, 2)==0)
        pa(i*2+1:i*2+2) = [0; d];
        pb(i*2+1:i*2+2) = [0; -d];
    else
        pa(i*2+1:i*2+2) = [0; -d];
        pb(i*2+1:i*2+2) = [0; d];
    end 
end 

% make bodies
bodies = [];
constraints = [];
parent = -1;
z = {};
for i = 0:numBod-1
    Mi = M(i*3+1:i*3+3, i*3+1:i*3+3);
    if i == numBod -1
        children = [];
        b = make_test_body(i+1, Mi, parent, children);
    else
        children = [numBod + i + 1];
        c = make_constraint(numBod + i + 1, zeros(2, 3), i+1, i+2);
        constraints = [constraints c];
        b = make_test_body(i+1, Mi, parent, children);
        parent = children(1);
    end
    bodies = [bodies b];
    z = [z; zeros(3, 1)];
end

allnodes = [bodies constraints];
for i=1:size(constraints, 2)-1
    z = [z; zeros(2, 1)];
end

%% Run 
q = q0; % initial position of the system
qdot = zeros(numBod*3,1); % initial velocity of the system
%           w vx vy ...
allConstrF = zeros(numBod*3,stps);
% external force only applies for an initial period of time
extF = zeros(numBod*3,stps);
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
    J = Jfuncmod(pa, pb, q);

    [allCOM, allBars, allax] = getAllBars(q,w,h,j,d,m);
    visualizeAllBars(allCOM, allBars, allax);

    c = cfuncmod(pa, pb, q, qdot);

    b = -(J*inv(M)*extF(:, i) + c);
    
    for bi=0:size(constraints, 2)-1
        allnodes(numBod+bi+1).D = J(bi*2+1:bi*2+2, bi*3+1:bi*3+6);
        z{numBod+bi+1} = -b(bi*2+1:bi*2+2);
    end

    % A\b NAIEVE SOLVE
    A = Afunc(pa, pb, q);
    lambda = A\b; %inv(A)*b;

    % SPARSE SOLVE
%     [H, forwards] = sparsefactor(allnodes);
%     ylamb = sparsesolve(H, z, allnodes, forwards);
%     lambda = cell2mat(ylamb(size(bodies, 2)+1:end));

    % DENSE SOLVE
%     [H, H_tree, forwards] = densefactor(allnodes);
%     ylamb = densesolve(H, z, forwards);
%     xdense = cell2mat(ylamb);
%     lambda = cell2mat(ylamb(size(bodies, 2)+1:end));
    
    qddot = inv(M)*J'*lambda + inv(M)*extF(:, i); % (6x1) update the velocity
    allConstrF(:,i) = J'*lambda;

    [q, qdot] = forwardeuler(q, qdot, qddot, dt);

    frame = getframe(gcf);
    writeVideo(vid,frame);
    cla
end
close(vid);