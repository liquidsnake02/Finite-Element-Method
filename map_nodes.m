%% maps the contributions from each element to the appropriate locations
%  in the matrix, so that the contribution can be added to the total matrix
%  input: 
%  the dimensions of the total matrix, 
%  the element whose contribution is being added
%  the local contribution of the element
%  output: 
%  the contribution to the total matrix
function contribution = map_nodes(dimensions, element, local, is_force)
    nodes = element.nodes;
    contribution = zeros(dimensions);
    if is_force == 1
        contribution = zeros(dimensions, 1);
        for j = 1:length(local)
            contribution(nodes(j)) = local(j);
        end
    else
        for j = 1:length(local)
            for k = 1:length(local)
                contribution(nodes(j),nodes(k)) = local(j,k);
            end
        end
    end
end