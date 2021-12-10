clear all;

% two bar example
% define all constriant blocks and mass blocks
% define tree structure too

% define N bodies
body1.i = 1;
body1.dim = 3;
body1.isConstraint=0;
body1.D = eye(3);
body1.parent = -1; % on parent
body1.children = [3];

body2.i = 2;
body2.dim = 3;
body2.isConstraint=0;
body2.D = eye(3) * 2; % mass matrix
body2.parent = 3;
body2.children = [];

% define M constraints
c1.i = 3;
c1.dim = 2;
c1.isConstraint=1;
c1.D = [1, 0, 1, 2, 0, 2;
        0, -1, 0, 0, -2, 0]; % jacobian
c1.parent = 1;
c1.children = [2];

bodies = [body1, body2];
constraints = [c1];
allnodes = [bodies constraints];
y = [2;3];
b = {zeros(3, 1); zeros(3, 1); y};

% direct solve
M = zeros(6, 6);
M(1:3, 1:3) = body1.D;
M(4:6, 4:6) = body2.D;
J = zeros(2, 6);
J = c1.D;

H_d = zeros(8, 8);
H_d(1:6, 1:6) = M;
H_d(7:8, 1:6) = J;
H_d(1:6, 7:8) = J';

b_d = cell2mat(b);
x_d = H_d \ b_d

% dense solve
[H, H_tree, forwards] = densefactor(bodies, constraints, allnodes);
x = densesolve(H, b, forwards);
cell2mat(x)

% sparse solve
[H, forwards] = sparsefactor(bodies, constraints, allnodes);
x = sparsesolve(H, b, allnodes, forwards);
cell2mat(x)
