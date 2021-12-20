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

%% System definition! Change all parameters here.

hTree = 10; % height of the binary tree
% 2^0 + 2^1 + ... + 2^h = 2^(h+1)-1
totBod = 2^(hTree+1)-1;

root = [-pi/2,0,0]';
relth =  zeros(totBod-1,1); % initial position of the system (6x1)
%th = pi/4;
q0 = zeros(totBod*3,1);%generateInitCoords(root, relth, d);
%q0 = generateInitTreeCoords(root, th, d, hTree);
numBod = totBod; length(q0)/3;

% Apply a tip force
Fext = zeros(numBod*3,1);
%Fext(end-2:end) = [0 -25 50]';

dt = 0.025; % timestep size [s]
time = 1; % total simulation time [s]

solverType = 1; % 1: A\b, 2: sparse, 3: dense
auxConstraint = 1; % 1: no aux constraints, 2: auxillary constraints
visualize = 0; % 1=visualize, otherwise=no

%% Get the Jacobian

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

for i = 0:hTree -1 % for every level of the tree
    for k = 1:2^i % for every block node
        ind = (k-1)*2;
        par = 2^i+k-1; % parent

        ln = 2^(i+1)+ind;
        rn = 2^(i+1)+ind+1;

        % left branch            
        xa = cat(1,xa,R2(q(par,1))*pa(par,:)' + q(par,2:3)');
        xb = cat(1,xb,R2(q(ln,1))*pb(ln,:)' + q(ln,2:3)');

        % right branch
        xa = cat(1,xa,R2(q(par,1))*pa(par,:)' + q(par,2:3)');
        xb = cat(1,xb,R2(q(rn,1))*pb(1,:)' + q(rn,2:3)');
    end
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

%% Jacobian of aux constraint

% constrain the COM of the first joint
aend = R2(q(1)) * [pa(1); -pa(2)] + q(2:3);
a = root(2:3) + R2(root(1))*[0; -d];

Ca = aend - a;
Ja = jacobian(Ca, q);

Ca_dot = Ja * qdot;
ca_sym = jacobian(Ca_dot, q) * qdot + ks * Ca + kd * Ca_dot;
Ca_ddot = Ja * qddot + jacobian(Ca_dot, q) * qdot;

Jafunc = matlabFunction(Ja, 'Vars', {pa, pb, q, qdot});
ca_symfunc = matlabFunction(ca_sym, 'Vars', {pa, pb, q, qdot});

%% Time Integration

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
%%
% make bodies
bodies = [];
constraints = [];
currparent = -1;
z = {};

% for tree structures

for i = 0:hTree % for every level of the tree
    Mi = M(i*3+1:i*3+3, i*3+1:i*3+3);
    
    if i==hTree % we are all the leaf of the tree
        children = [];
        for k = 1:2^(i)
            ind = 2^i+k-1;
            b = make_test_body(ind, Mi, numBod + ind -1, children);
            bodies = [bodies b];
            z = [z; zeros(3, 1)];
        end 

    else % another level of the tree
        
        for k = 1:2^i % for every block node
            ind = 2^i+k-1;
            children = [totBod + ind*2-1, totBod + ind*2];% there are two child constraints

            cl = make_constraint(totBod + ind*2-1, zeros(2, 3), ind, ind*2);
            cr = make_constraint(totBod + ind*2, zeros(2, 3), ind, ind*2+1);

            constraints = [constraints cl cr];
           
            currparent = numBod+ind-1;
            if(i==0)
                currparent = - 1;
            end
            b = make_test_body(ind, Mi, currparent, children);
            bodies = [bodies b];
            z = [z; zeros(3, 1)];
        end 
    end 
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

if visualize == 1
    vid = VideoWriter('video.avi');
    open(vid);
    figure(1)
    hold on
    xlim([-15 15])
    ylim([-15 15])
    grid on
end

tStart = tic;
for i=1:size(tall, 2)

    i

    t = tall(i);

    [allCOM, allBars, allax] = getAllBars(q,w,h,j,d,m);
    if visualize == 1
        visualizeAllBars(allCOM, allBars, allax);
    end

    % Get J, c, b for primary constraint
    J = Jfuncmod(pa, pb, q);
    c = cfuncmod(pa, pb, q, qdot);
    b = -(J*inv(M)*extF(:, i) + c);
    
    % (1) solve lambda for primary constraint
    if(solverType == 1) % (1a) A\b NAIEVE SOLVE
        tic
        A = Afunc(pa, pb, q);
        lambda = A\b; %inv(A)*b;
        toc
    elseif(solverType == 2) % (1b) SPARSE SOLVE
        for bi=0:size(constraints, 2)-1
            par = floor(bi/2);
            child = bi+1;
            Jic = J(bi*2+1:bi*2+2, par*3+1:par*3+3);
            Jip = J(bi*2+1:bi*2+2, child*3+1:child*3+3);
            allnodes(numBod+bi+1).D = [Jic Jip];
            z{numBod+bi+1} = -b(bi*2+1:bi*2+2);
        end
        tic
        [H, forwards] = sparsefactor(allnodes);
        ylamb = sparsesolve(H, z, allnodes, forwards);
        toc
        lambda = cell2mat(ylamb(size(bodies, 2)+1:end));

    elseif(solverType == 3) % (1c) DENSE SOLVE
        for bi=0:size(constraints, 2)-1
            par = floor(bi/2);
            child = bi+1;
            Jic = J(bi*2+1:bi*2+2, par*3+1:par*3+3);
            Jip = J(bi*2+1:bi*2+2, child*3+1:child*3+3);
            allnodes(numBod+bi+1).D = [Jic Jip];
            z{numBod+bi+1} = -b(bi*2+1:bi*2+2);
        end
        
        tic
        [H, H_tree, forwards] = densefactor(allnodes);
        ylamb = densesolve(H, z, forwards);
        toc
        xdense = cell2mat(ylamb);
        lambda = cell2mat(ylamb(size(bodies, 2)+1:end));
    end 
    
    if(auxConstraint == 1)
        % no aux constraints
        qddot = inv(M)*J'*lambda + inv(M)*extF(:, i); % (6x1) update the velocity
    elseif(auxConstraint == 2)
        % (2) evaluate Ja K ca
        vdot_aux = inv(M) * (J'*lambda + extF(:, i));
        Ja = Jafunc(pa, pb, q, qdot);
        K = Ja';
        ca = ca_symfunc(pa, pb, q, qdot);

        MhatK = zeros(size(K));

        if(solverType == 1)
            % (3) solve for M_hat * K
            for k=1:size(K, 2)
                b = -J * inv(M) * K(:, k);
                lambda = A\b;
                MhatK(:, k) = inv(M) * (J' * lambda + K(:, k));
            end
        
            % (4) solve for mu (O (K^3))
            mu = (Ja*MhatK) \ (a - Ja*vdot_aux - ca);
        
            % (5) solve for primary response to Fext + K*mu
            b = -(J*inv(M)*(K*mu + extF(:, i)) + c);
            lambda = A\b;
            
            % aux constraints
            qddot = inv(M)*(K*mu + J'*lambda + extF(:, i)); % (6x1) update the velocity
        else % Sparse or Dense Solve
            % (3) solve for M_hat * K
            MhatK = zeros(size(K));
            for k=1:size(K, 2)
                b = -J * inv(M) * K(:, k);
        
                for bi=0:size(constraints, 2)-1
                    z{numBod+bi+1} = -b(bi*2+1:bi*2+2);
                end

                ylamb = sparsesolve(H, z, allnodes, forwards);
                lambda = cell2mat(ylamb(size(bodies, 2)+1:end));
                MhatK(:, k) = inv(M) * (J' * lambda + K(:, k));
            end
        
            % (4) solve for mu (O (K^3))
            mu = (Ja*MhatK) \ (a - Ja*vdot_aux - ca);
            
            % (5) solve for primary response to Fext + K*mu
            b = -(J*inv(M)*(K*mu + extF(:, i)) + c);
            for bi=0:size(constraints, 2)-1
                z{numBod+bi+1} = -b(bi*2+1:bi*2+2);
            end

            ylamb = sparsesolve(H, z, allnodes, forwards);
            lambda = cell2mat(ylamb(size(bodies, 2)+1:end));
        end

        % aux constraints
        qddot = inv(M)*(K*mu + J'*lambda + extF(:, i)); % (6x1) update the velocity
    end

    allConstrF(:,i) = J'*lambda;

    [q, qdot] = forwardeuler(q, qdot, qddot, dt);

    if visualize == 1
        frame = getframe(gcf);
        writeVideo(vid,frame);
        cla
    end
end
elapsed = toc(tStart);
elapsed = elapsed / size(tall, 2)
if visualize == 1
    close(vid);
end