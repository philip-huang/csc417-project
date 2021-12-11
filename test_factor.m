%% Test case 1
clear all;

% two bar example
% define all constriant blocks and mass blocks
% define tree structure too

% define N bodies
body1.i = 1;
body1.dim = 3;
body1.isConstraint=0;
body1.D = eye(3);
body1.D(2, 2) = 0.1;
body1.parent = -1; % on parent
body1.children = [3];

body2.i = 2;
body2.dim = 3;
body2.isConstraint=0;
body2.D = eye(3) * 2; % mass matrix
body2.D(2, 2) = 0.3;
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
y = [-2;10];
b = {zeros(3, 1); zeros(3, 1); y};

% direct solve
M = zeros(6, 6);
M(1:3, 1:3) = body1.D;
M(4:6, 4:6) = body2.D;
J = zeros(2, 6);
J = c1.D;

H_d = zeros(8, 8);
H_d(1:6, 1:6) = M;
H_d(7:8, 1:6) = -J;
H_d(1:6, 7:8) = -J';

b_d = cell2mat(b);
x_d = H_d \ b_d

% dense solve
[H, H_tree, forwards] = densefactor(allnodes);
x = densesolve(H, b, forwards);
cell2mat(x)

% sparse solve
[H, forwards] = sparsefactor(allnodes);
x = sparsesolve(H, b, allnodes, forwards);
cell2mat(x)

%% Test case 2
clear all;
% 6 bar example
b1 = make_test_body(1, eye(3), -1, [7, 10]);
b2 = make_test_body(2, eye(3)*2, 7, [8, 9]);
b3 = make_test_body(3, eye(3)*1.3, 8, []);
b4 = make_test_body(4, eye(3)*1.4, 9, []);
b5 = make_test_body(5, eye(3)*3, 10, [11]);
b6 = make_test_body(6, eye(3)*2.2, 11, []);
J = [1, 0, 1 2, 0, -2;
    0   -1 0 0  3   0];
c1 = make_constraint(7, J, 1, [2]);
c2 = make_constraint(8, J, 2, [3]);
c3 = make_constraint(9, J, 2, [4]);
c4 = make_constraint(10, J, 1, [5]);
c5 = make_constraint(11, J, 5, [6]);

bodies = [b1 b2 b3 b4 b5 b6];
constraints = [c1 c2 c3 c4 c5];
allnodes = [bodies constraints];
b = {};
for i=1:numel(bodies)
    b{i, 1} = zeros(3, 1);
end
b{7, 1} = [10; 1];
b{8, 1} = [2; 3];
b{9, 1} = [0; 0];
b{10, 1} = [11; 1];
b{11, 1} = [-4; -1];

% direct solve
[M, J] = get_M_J(bodies, constraints);
H = [M, -J';
    -J, zeros(size(J, 1))];
bd = cell2mat(b);
x_direct = H \bd;

% dense solve
[H, H_tree, forwards] = densefactor(allnodes);
x = densesolve(H, b, forwards);
xdense = cell2mat(x);

% sparse solve
[H, forwards] = sparsefactor(allnodes);
x = sparsesolve(H, b, allnodes, forwards);
xsparse = cell2mat(x);

assert(all(xdense==xsparse))

function [body] = make_test_body(i, M, p, c)
    body.i = i;
    body.dim = size(M, 1);
    body.isConstraint = 0;
    body.D = M;
    body.parent = p;
    body.children = c;
end

function [constraint] = make_constraint(i, J, p, c)
    constraint.i = i;
    constraint.dim = size(J, 1);
    constraint.isConstraint = 1;
    constraint.D = J;
    constraint.parent = p;
    constraint.children = c;
end
