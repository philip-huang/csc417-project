function [M,J] = get_M_J(bodies, constraints)
% get the big matrix M and J
M = {};
for i=1:numel(bodies)
    for j=1:numel(bodies)
        if i==j
            M{i, i} = bodies(i).D;
        else
            M{i, j} = zeros(bodies(i).dim, bodies(j).dim);
        end
    end
end

M=cell2mat(M);

J = {};
for i=1:numel(constraints)
    p = constraints(i).parent;
    c = constraints(i).children(1);
    Ji = constraints(i).D;
    dimi = constraints(i).dim;
    dimp = bodies(p).dim;
    for j=1:numel(bodies)
        if j==p
            J{i, j} = Ji(:, 1:dimp);
        elseif j==c
            J{i, j} = Ji(:, dimp+1:end);
        else
            J{i, j} = zeros(dimi, bodies(j).dim);
        end
    end
end
J = cell2mat(J);

end

