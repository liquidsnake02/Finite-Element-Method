%% performs gauss quadrature on an arbitrary function
%  input: 
%  an anonymous function F
%  gauss points gp
%  gauss weights gw
%  output: 
%  the approximate integral of F 
function q = quadrature(F, gp, gw)
q = 0;
    for j=1:length(gp)
        q = q + gw(j)*F(gp(1,j), gp(2,j));
    end
end