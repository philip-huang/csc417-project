function [x] = sparsesolve(H, z, allnodes, forwards)
%DENSESOLVE

n = numel(forwards);
% shift z
xx = cell(size(z));
for i=1:numel(forwards)
    xx{i} = z{forwards(i)};
end

% solve
for i=1:n
    node = allnodes(i);
    for jind = 1:numel(node.children)
        j = node.children(jind);
        xx{i} = xx{i} - H{i, j} *xx{j};
    end
end
for i=n:-1:1
    xx{i} = inv(H{i, i}) * xx{i};
    if i ~= n
        p = allnodes(i).parent;
        xx{i} = xx{i} - H{i, p} * xx{p};
    end
end

% shift x back
x = cell(size(xx));
for i=1:numel(forwards)
    x(forwards(i)) = xx(i);
end


