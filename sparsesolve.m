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
    node = allnodes(forwards(i));
    for jind = 1:numel(node.children)
        j = find(forwards==node.children(jind));
        xx{i} = xx{i} - H{j, i}' *xx{j};
    end
end
for i=n:-1:1
    node = allnodes(forwards(i));
    xx{i} = H{i, i} \ xx{i};
    if i ~= n
        p = find(forwards==node.parent);
        xx{i} = xx{i} - H{i, p} * xx{p};
    end
end

% shift x back
x = cell(size(xx));
for i=1:numel(forwards)
    x(forwards(i)) = xx(i);
end


