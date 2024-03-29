function [H, H_tree, forwards] = densefactor(allnodes)
% Get H eq (10) from M and J
% Input M: mass matrix
% Input J: Jacobian matrix
% dofb: dof of the body
% dofc: dof of the constraint
% return H

[forwards, backwards] = dfs(allnodes(1), allnodes, [], []);

% dim of block
H = {};

% mass block first
for i=1:numel(forwards)
    node = allnodes(forwards(i));
    for j=1:numel(forwards)
        if node.isConstraint == 0
            if i==j
                H{i, j} = node.D;
            else
                H{i, j} = zeros(node.dim, allnodes(forwards(j)).dim);
            end
        end
    end
end

for i=1:numel(forwards)
    node = allnodes(forwards(i));
    for j=1:numel(forwards)
        if node.isConstraint == 1
            if node.parent == forwards(j)
                dimp = allnodes(node.parent).dim;
                Jip = node.D(:, 1:dimp);
                H{i, j} = -Jip;
                H{j, i} = -Jip';
            elseif node.children(1) == forwards(j)
                dimp = allnodes(node.parent).dim;
                Jic = node.D(:, dimp+1:end);
                H{i, j}= -Jic;
                H{j, i} = -Jic';
            else
                H{i, j} = zeros(node.dim, allnodes(forwards(j)).dim);
            end
        end
    end
end
H_tree = H;

%cell2mat(H)   

% densefactor
for i=1:numel(forwards)
    for k=i-1:-1:1
        H{i, i} = H{i, i} - H{k, i}'*H{k, k}*H{k, i};
    end
    for j=i+1:numel(forwards)
        for k=i-1:-1:1
            H{i, j} = H{i, j} - H{k, i}'*H{k, k}*H{k, j};
        end
        H{i, j} = inv(H{i, i})*H{i, j};
    end
end
%cell2mat(H)  

% debug
% Hf = cell2mat(H);
% D = diag(diag(Hf));
% Lt = triu(Hf, 1) + eye(size(Hf));
% L = tril(Hf, -1) + eye(size(Hf));
% Lt' * D * Lt - cell2mat(H_tree);
end

function [forwards, backwards] = dfs(node, allnodes, forwards, backwards)
    children = node.children;
    for c=1:numel(children)
        [forwards, backwards] = dfs(allnodes(children(c)), allnodes, forwards, backwards);
    end
    forwards = [forwards, node.i];
    backwards = [node.i backwards];
end