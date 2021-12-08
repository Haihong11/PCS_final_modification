function g=piecewise_expmap(x,theta,xi) %Exponential map

xihat          = vector_hat(xi);
if theta==0
    g           = eye(4)+x*xihat;
    
else
    g           = eye(4)+x*xihat+...
                 ((1-cos(x*theta))/(theta^2))*xihat^2+...
                 ((x*theta-sin(x*theta))/(theta^3))*xihat^3;
end