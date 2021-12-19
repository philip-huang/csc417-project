function [body] = make_test_body(i, M, p, c)
    body.i = i;
    body.dim = size(M, 1);
    body.isConstraint = 0;
    body.D = M;
    body.parent = p;
    body.children = c;
end