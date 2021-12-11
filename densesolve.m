function [x] = densesolve(H, z, forwards)
%DENSESOLVE

n = numel(forwards);
% shift z
xx = cell(size(z));
for i=1:numel(forwards)
    xx{i} = z{forwards(i)};
end

% solve
for i=1:n
    for j = 1:i-1
        xx{i} = xx{i} - H{j, i}' *xx{j};
    end
end
for i=n:-1:1
    xx{i} = inv(H{i, i}) * xx{i};
    for j=i+1:n
        xx{i} = xx{i} - H{i, j} * xx{j};
    end
end

% shift x back
x = cell(size(xx));
for i=1:numel(forwards)
    x(forwards(i)) = xx(i);
end


