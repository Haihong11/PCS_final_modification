function matrix= vector_hat(screw)

matrix         = zeros(4,4);
matrix(1:3,1:3)= vector_tilde(screw(1:3)); 
matrix(1:3,4)  = screw(4:6);