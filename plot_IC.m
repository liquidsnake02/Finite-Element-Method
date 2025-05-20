%% A script to plot the initial condition. 
%  Used for creating the plot in the paper
xd = [0:0.5:2*pi];
yd = [0:0.5:2*pi];
x0 = 2*pi;
y0 = x0;
[X,Y] = meshgrid(xd,yd);
g = @(x,y) 5.*sin(x/2 + pi).*cos((y + pi)/2);
G = g(X,Y);
surf(X,Y,G)
axis([0 2*pi 0 2*pi 0 5])
title('Plotting the intial condition on [0,2pi]x[0,2pi]')
xlabel('x')
ylabel('y')
zlabel('g(x,y)')
figure
contour(X,Y,G)
axis([0 2*pi 0 2*pi 0 5])
title('Plotting the intial condition on [0,2pi]x[0,2pi]')
xlabel('x')
ylabel('y')
zlabel('g(x,y)')