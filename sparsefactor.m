function [H, forwards] = sparsefactor(bodies,constraints, allnodes)
% Get H eq (10) from M and J
% Input M: mass matrix
% Input J: Jacobian matrix
% return H

% allnodes1 should be the root node
[forwards, backwards] = dfs(allnodes(1), allnodes, [], []);
n = numel(forwards);

% dim of block
H = {};

for i=1:numel(forwards)
    node = allnodes(forwards(i));
    if node.isConstraint
        p = find(forwards==node.parent);
        c = find(forwards==node.children(1));
        
        dimp = allnodes(node.parent).dim;
        Jip = node.D(:, 1:dimp);
        Jic = node.D(:, dimp+1:end);
        H{i, p} = Jip;
        H{p, i} = Jip';
        H{i, c} = Jic;
        H{c, i} = Jic';
        H{i, i} = zeros(node.dim);
    else
        H{i, i} = node.D;
    end
end

    for i = 1:numel(forwards)
    node = allnodes(forwards(i));
    for kind = 1:numel(node.children)
        k = find(forwards==node.children(kind));
        H{i, i} = H{i, i} - H{k, i}'*H{k, k} * H{k, i};
    end
    if i ~= n
        p = find(forwards==node.parent);
        H{i, p} = inv(H{i, i}) * H{i, p};
    end
end
    
end

function [forwards, backwards] = dfs(n, allnodes, forwards, backwards)
    children = n.children;
    for c=1:numel(children)
        [forwards, backwards] = dfs(allnodes(children(c)), allnodes, forwards, backwards);
    end
    forwards = [forwards, n.i];
    backwards = [n.i backwards];
end