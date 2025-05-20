%% generates a mesh over a square domain
%  input: h is the stepsize that controlls how fine the mesh  
%         For the following descriptions, view the plane with the origin
%         in the bottom left
%         x1 is the x-coordinate of the bottom left corner
%         x2 is the x-coordinate of the bottom right corner
%  output: an array of triangle structs that covers the square mesh and the
%          set of points that lie on the boundary
function [tri_array, boundary_set] = generate_mesh(h, x1, x2)
    % sqare domain, so we know the height from the width
    % cut the domain into pieces
    xsides = [x1:h:x2];
    % reassemble the square matrix of ordered pairs
    pointsperside = length(xsides);
    zeropoints = zeros(1,length(xsides));
    xpoints = [xsides ; zeropoints];
    
%% generates the mesh as a struct array of triangle objects
trianglesperside = pointsperside - 1;
% creates empty triangle objects and adds it to the last cell of the array
% as a way to preallocate array size
triangle = struct('points', [0 0 0; 0 0 0],...
    'nodes', [0,0,0], 'boundary', [0 0 0]);
triangles(trianglesperside, 2*trianglesperside) = triangle;
tricol(trianglesperside,1) = triangle;
triangles(:,1) = tricol;

% loops through the triangles adding them to the array 
tricol_a(trianglesperside,1) = triangle;
tricol_b(trianglesperside,1) = triangle;
boundary_set = zeros(4*pointsperside - 4, 1);
bindex = 1;
for n = 1:trianglesperside
    % build tricol_a and tricol_b
    abase1 = xpoints(:,n);
    abase2 = xpoints(:,n+1);
    bbase1 = abase1 + [0 h]';
    bbase2 = abase2 + [0 h]';
    for m = 1:trianglesperside
        nodes = [(m-1)*pointsperside+n,...
        (m-1)*pointsperside+n+1, m*pointsperside + n];
        [boundary, boundary_set, bindex] =...
            check_boundary(nodes, pointsperside, boundary_set, bindex);
        tricol_a(m) = struct('points', [abase1, abase2, bbase1],...
        'nodes', nodes, 'boundary', boundary);

        nodes = [m*pointsperside + n + 1,...
        m*pointsperside + n, (m-1)*pointsperside+n+1];
        [boundary, boundary_set, bindex] =...
            check_boundary(nodes, pointsperside, boundary_set, bindex);
        tricol_b(m) = struct('points', [bbase2, bbase1, abase2],...
        'nodes', nodes, 'boundary', boundary);
   
        abase1 = bbase1;
        abase2 = bbase2;
        bbase1 = bbase1 + [0 h]';
        bbase2 = bbase2 + [0 h]';
    end
    triangles(:,2*n-1) = tricol_a;
    triangles(:,2*n) = tricol_b;
end
% return the generated mesh
tri_array = triangles;
end
function [boundary, boundary_set, bindex]=...
    check_boundary(nodes, pointsperside, boundary_set, bindex)
    boundary = [0 0 0];
    i = 1;
    for n = nodes
        if (n < pointsperside)
            boundary(i) = 1;
            [boundary_set, bindex] =...
                add_boundary_node(n, boundary_set, bindex);
        end
        if (mod(n, pointsperside) == 0 || mod(n, pointsperside) == 1)
            boundary(i) = 1;
            [boundary_set, bindex] =...
                add_boundary_node(n, boundary_set, bindex);
        end
        if (n > pointsperside*(pointsperside - 1))
            boundary(i) = 1;
            [boundary_set, bindex] =...
                add_boundary_node(n, boundary_set, bindex);
        end
        i = i + 1;
    end
end
function [boundary_set, boundary_index] =...
    add_boundary_node(node, bset, bindex)
    for j = 1:bindex
        if (bset(j) == node)
            boundary_set = bset;
            boundary_index = bindex;
            return;
        end
    end
    bset(bindex) = node;
    boundary_set = bset;
    boundary_index = bindex + 1;
end