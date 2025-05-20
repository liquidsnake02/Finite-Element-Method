%% solves the following heat equation with the Fourier series solution
%  0<= x,y <= 2*pi
%  u(x,0,t) = u(x,2pi,t) = 0
%  u(0,y,t) = u(2pi,y,t) = 0
%  u(x,y,0) = 5sin(x/2 + pi)cos((y+pi)/2)
%  input: 
%  run is a number that controls whether the function runs or not.
%  if run = 0 then the function does not run, otherwise it does.
%  tspan is the period of time that the function runs over. more than ten
%  discrete points in time leads to incredibly long run time
%  plot_contour decides whether contour plots are created. They will plot
%  for all values except 0
function fourier_solution(run, sum_limit, tspan, dx, plot_contour)
tic % measuring runtime
if(run == 0)
    return;
end
T=tspan;
N = sum_limit;
M = sum_limit;
xd = [0:dx:2*pi];
yd = [0:dx:2*pi];
x0 = 2*pi;
y0 = x0;
pcontour = plot_contour;
[X,Y] = meshgrid(xd,yd);
f = @(x,y) 5.*sin(x/2 + pi).*cos((y + pi)/2);
v = 4/(x0*y0);
u = zeros(length(X));
for t = T
    for m = 1:M
        for n = 1:N
            intf = @(x,y) f(x,y).*sin(m.*pi.*x./x0).*sin(n.*pi.*y./y0);
            Amn = v*integral2(intf,0,x0,0,y0);
            e = exp(( -pi^2 * (m^2/(x0^2) + n^2/(y0^2)) ).*t);
            u = u + Amn.*sin(m.*pi.*X./x0).*sin(n.*pi.*Y./y0).*e;
        end
    end
    figure
    surf(X,Y,u)
    axis([0 2*pi 0 2*pi 0 5])
    title("The solution to the equation at time t = " + t)
    xlabel('x')
    ylabel('y')
    zlabel("u(x,y," + t + ")")
    if pcontour ~= 0
        figure
        contour(X,Y,u)
        axis([0 2*pi 0 2*pi 0 5])
        title("The solution to the equation at time t = " + t)
        xlabel('x')
        ylabel('y')
        zlabel("u(x,y," + t + ")")
    end
    u = reshape(u, [], 1);
    u = zeros(length(X));
end
fprintf("Elapsed Fourier time is %d\n", toc)
end