%% The main file for running the FEM implementation and the Fourier
%  solution.

%% Solution selection
%  The FEM and Fourier variables will decide if the respective solution
%  runs. The solutions will run for any value that is not 0, so you must
%  set them to 0 to restrict them. The Fourier solution is very slow, so I
%  recommend leaving it off unless you have a reason to look at it.
%
%  time_points decides which time steps are plotted by the program for both
%  the FEM and Fourier solutions.
%
%  dx is the spatial step for both solutions. Be wary of setting this
%  too high, as it can drastically increase runtime.
%
%  contour will decide if the solutions plot contour plots as well as
%  surface plots. Contours will plot for any value that is not 0. If you 
%  are doing running both solutions with multiple key time points then I 
%  recommend setting this to 0.
FEM = 1;
Fourier = 1;
time_points = [0 1];
dx = 0.5;
contour = 0;

%% FEM section
%  left and right bound decide the region of the problem. It will be a
%  square with those endpoints
%
%  initial_condition is exactly what it sounds like. It must be an
%  anonymous function handle, so if you decide to change it then do not
%  remove the @(x,y).
%
%  forces is the source term. Set to 0 for a homogeneous problem. This also
%  must be an anonymous function handle, so do not remove @(x,y)
%
%  dt is the timestep being used for implicit Euler's
%
%  remove_negtaives will set all negative values in a solution to 0 when
%  plotting it. Setting it to 0 does not remove the negatives, any other 
%  value will
left_bound = 0;
right_bound = 2*pi;
initial_condition = @(x,y) 5*sin(x/2 + pi)*cos((y+pi)/2);
forces = @(x,y) 0;
dt = 0.01;
remove_negatives = 1;

FEM_2D_HeatEq(FEM, left_bound, right_bound, dx,...
    initial_condition, forces, time_points, dt, remove_negatives, contour);

%% Fourier section
%  The fourier solution of the main problem that appears in section 6.1 of
%  the paper.
%
%  sum_limit is how high the double fourier series sum indices will go. the
%  fourier solution here is computed using matlab numerical integration,
%  and so a high sum limit will heavily increase run time. The base value
%  is 10, and that will still take a noticable amount of time to compute.
sum_limit = 10;

fourier_solution(Fourier, sum_limit, time_points,...
    dx, contour);