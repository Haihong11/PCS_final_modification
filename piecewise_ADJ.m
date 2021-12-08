%--------Tangent operator of the exponential map-----%

function ADJg=piecewise_ADJ(x,theta,xi)  

adjxi       =matrix_adj(xi);

if theta==0
    ADJg        =x*eye(6)+((x^2)/2)*adjxi;
    
else
    ADJg        =x*eye(6)+((4-4*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^2))*adjxi+...
                 ((4*x*theta-5*sin(x*theta)+x*theta*cos(x*theta))/(2*theta^3))*adjxi^2+...
                 ((2-2*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^4))*adjxi^3+...
                 ((2*x*theta-3*sin(x*theta)+x*theta*cos(x*theta))/(2*theta^5))*adjxi^4;
end