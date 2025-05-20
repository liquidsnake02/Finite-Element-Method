%% Solves the given 2D heat problem with the finite element method
%  input: 
%  run decides if the program runs. Will run for any nonzero value
%  left and right bound and the dimensions of the problem
%  dx is the spatial step size
%  initial condition is the intial condition as a function
%  source is the source term as a function
%  time_points is the points in time being plotted
%  dt is the temporal step size
%  remove_negatives will remove all negatives from a soltuion unless it is
%  set to 0
%  plot_contour will plot contour plots along with the surfaces unless it
%  is set to 0
function FEM_2D_HeatEq(run, left_bound, right_bound, dx,...
    initial_condition, source, time_points, dt, remove_negatives, plot_contour)
tic % measuring runtime
%%  Check run
if (run == 0)
    return;
end
%% set problem and boundaries
x0 = left_bound;
x1 = right_bound;
h = dx;
domain = [x0:h:x1];
domaindim = length(x0:h:x1);
numofnodes = domaindim^2;
IC = initial_condition;
forces = source;
tstep = dt;
T = time_points;
rz = remove_negatives;
pcontour = plot_contour;
%% generate spatial mesh and adjacent properties
%  creates a uniform mesh where each h x h square is split by a diagonal,
%  resulting in two right triangles per square
[mesh, boundary_set] = generate_mesh(h,x0,x1);
%  defining the shape function, which forms a triangle base pyramid over
%  the reference element. each piece of the shape function is 1 at only a
%  single "node" of the reference element, and 0 at all others.
%  phi = phi1 + phi2 + phi3
phi = {@(x,y) 1-x-y, @(x,y) x, @(x,y) y};
%  we will need a jacobian matrix to perform our coordinate transforms for
%  each element to map to the reference element. the jacobian can be
%  calculated as [differential matrix]*[point matrix]'. here we will enter
%  the differential matrix by hand because it is independent of the points
%  being tranformed
phi1dx = -1;
phi1dy = -1;
phi2dx = 1;
phi2dy = 0;
phi3dx = 0;
phi3dy = 1;
gradphi = [phi1dx phi2dx phi3dx; phi1dy phi2dy phi3dy];
%  jacobian function, also generates the inverse because both are needed
%  frequently
jacobian = @(triangle) gradphi*triangle.points';
%  we will use gauss quadrature to calculate integrals using three gauss
%  points. for a triangle, which we will hard code in
gauss_points = [1/6 3/2 1/6; 0 1/6 2/3];
gw = [1/6 1/6 1/6];
%% create linear system
%  same dimensions as domain bc 1 degree of freedom on all elements
%  stiffness matrix S 
stiffness = zeros(numofnodes); 
%  mass matrix M
mass = zeros(numofnodes);
%  force vector f
force = zeros(numofnodes,1);
for i = 1:length(mesh)/2
    for j = 1:length(mesh)
        J = jacobian(mesh(i,j));
         invJ = inv(J);
         local_stiffness = zeros(3); % 3 bc using triangles
         local_mass = zeros(3);
         local_force = zeros(3,1);
        for n = 1:3
            forceseq = @(x,y) forces(x,y)*phi{n}(x,y)*abs(det(J));
            local_force(n)= quadrature(forceseq, gauss_points, gw);
            for m = 1:3
            stiffnesseq =...
                @(x,y) gradphi(:,n)'*(invJ'*invJ)*gradphi(:,m)*abs(det(J));
            local_stiffness(n,m)=quadrature(stiffnesseq, gauss_points, gw);
            masseq =...
                @(x,y) phi{n}(x,y)*phi{m}(x,y)*abs(det(J));
            local_mass(n,m)=quadrature(masseq, gauss_points, gw);
            end
        end
        stiffness = stiffness +...
            map_nodes(numofnodes, mesh(i,j), local_stiffness, 0);
        mass = mass + map_nodes(numofnodes, mesh(i,j), local_mass, 0);
        force = force + map_nodes(numofnodes, mesh(i,j), local_force, 1);
    end
end
%  impose boundary on stiffness matrix
for b = boundary_set
    for d = 1:length(mass)
        stiffness(b,d) = 0;
    end
end
for b = 1:length(mass)
    for d = boundary_set
        stiffness(b,d) = 0;
    end
end
%% solve linear system
c_matrix = zeros(domaindim);
for i = 1:domaindim
    for j = 1:domaindim
        c_matrix(i,j) = IC(domain(i), domain(j));
    end
end
 %  plots the soltutions at t = 0
if (T(1) == 0)
    [X,Y] = meshgrid(domain,domain);
    figure
    surf(X,Y,c_matrix)
    axis([x0 x1 x0 x1 0 5])
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y,0)')
    title("FEM Approximation at time t = 0")
    if pcontour ~= 0
        figure
        contour(X,Y,c_matrix)
        axis([x0 x1 x0 x1 0 5])
        xlabel('x')
        ylabel('y')
        zlabel('u(x,y,0)')
        title("FEM Approximation at time t = 0")
    end
end
c_vector = reshape(c_matrix, [], 1);
for time_point = 1:length(T) - 1
    tspan = [T(time_point):tstep:T(time_point + 1)];
    %  implicit Euler up to time t
    for j = tspan
            rhs_vec = tstep*force + mass*c_vector;
            c_vector = (mass+tstep*stiffness)\rhs_vec;
    end
    %  reimpose lost BCs on solution
    for b = boundary_set
            c_vector(b) = 0;
    end
    if rz ~= 0
        for j = 1:length(c_vector)
            if c_vector(j) < 0
                c_vector(j) = 0;
            end
        end
    end
    c_matrix = reshape(c_vector, domaindim, domaindim);
    figure
    surf(X,Y,c_matrix)
    axis([x0 x1 x0 x1 0 5])
    xlabel('x')
    ylabel('y')
    zlabel("u(x,y," + tspan(end) + ")")
    title("FEM Approximation at time t = " + tspan(end))
    if pcontour ~= 0
        figure
        contour(X,Y,c_matrix)
        axis([x0 x1 x0 x1 0 5])
        xlabel('x')
        ylabel('y')
        zlabel("u(x,y," + tspan(end) + ")")
        title("FEM Approximation at time t = " + tspan(end))
    end
end
fprintf("Elapsed FEM time is %d\n", toc)
end