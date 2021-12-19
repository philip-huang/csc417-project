function [constraint] = make_constraint(i, J, p, c)
    constraint.i = i;
    constraint.dim = size(J, 1);
    constraint.isConstraint = 1;
    constraint.D = J;
    constraint.parent = p;
    constraint.children = c;
end