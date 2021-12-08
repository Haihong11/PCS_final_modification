%---------------------------Video--------------------------%
clear
clc
L         = [0.09,0.07,0.04];           % [m] one section of length of soft robot
num_piece = 3;                              % Number of piece [9,4,4]
num_disc  = [10,8,5];                        % Number of disc.
X1        = linspace(0,L(1),num_disc(1));   % coordinate of each disc;
X2        = X1(end) + linspace(0,L(2),num_disc(2));
X3        = X2(end) + linspace(0,L(3),num_disc(3));
X1(end)   = [ ];
X2(end)   = [ ];
X         = [X1,X2,X3];
dX        = L./(num_disc-1);                % The spacing between the discs
Rmax      = 0.01;                         % [m] Radius of base
Rmin      = 0.005;                         % [m] Radius of end

R_d       = (Rmax-Rmin) / sum(L);           % Radius reduction of the circular section along the X axis.

nsol      = 101;

%-------------get the solution--------------%

g_r        = [1 0 0 0; 
              0  1 0 0;
              0  0 1 0;
              0  0 0 1];

g_pre      = eye(4);


load('1.mat');

%----------Get the configuration of each disc------------%

for i = 1:nsol     
  
    Xi         = z(i,1:6*num_piece);
 
    g_pre      = eye(4);
    
    for ii = 1:num_piece
    
        xin          = Xi(6*ii-5:6*ii);                  
        kn           = xin(1:3);
        thetan       = sqrt(kn * kn');
        
        if ii == 1

            disc_index = 1;

        else

            disc_index = sum(num_disc(1:ii-1)) + 1;

        end
     
        num_disc_i   = num_disc(ii);
  
        for j = 1:num_disc_i 

            g(4*i-3:4*i,4*(disc_index-1)+4 * j-3:4*(disc_index-1) + 4 * j) ...
                                   = g_r * g_pre * piecewise_expmap(X(j),thetan,xin);
                               
        end

        g_pre        =  g_pre * piecewise_expmap(X(num_disc_i),thetan,xin);

    end

end

video         = VideoWriter(strcat('D:\PLS Cosserat model for multi-section soft manipulator dynamics\PCS_dynamic_ode45'));

FrameRate     = 10^2;      

open(video)

scrsz             = get(0,'ScreenSize'); 

figure('color','w','Position',[scrsz(3)/24 2*scrsz(4)/48 11*scrsz(3)/12 9*scrsz(4)/10])

axis equal
grid on
hold on

xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

%     view(0,0);
%     view(90,0);
%     view(45,45);
%     view(0,90);  
view(135,45);
% view(45,135);

axis ([-1.5*sum(L) 1.5*sum(L) -2*sum(L) 2*sum(L) -1.5*sum(L) 1.5*sum(L)])                   % XYZ in the Inertial frame  
% axis ([-num_piece*L num_piece*L -1.5*num_piece*L 1.5*num_piece*L -1.5*L 1.5*L])      % XYZ in the Inertial frame  

for i = 1:nsol

    cla      

    g1        = g(4*i-3:4*i,:);

    alpha     = linspace(0,2*pi,50);     % 100

    R_ii      = Rmax;

    for zz = 1:num_piece

        R_pre   = R_ii;
        
        if zz == 1

            disc_index = 1;

        else

            disc_index = sum(num_disc(1:zz-1)) + 1;

        end
         
        num_disc_i   = num_disc(zz);

        for ii = 1:num_disc_i

        %-------draw circle by the matrix, e.g., variable 'R'-------%  

            R_ii      = R_pre - R_d * (ii-1) * dX(zz);

            circle    = [zeros(1,50)         0;
                        0+R_ii * sin(alpha)  0;
                        0+R_ii * cos(alpha)  0;
                        ones(1,50)           1];     % (0,0) is origin of the circle.

            robot  = g1(:,4*(disc_index-1)+4*ii-3:4*(disc_index-1)+4*ii) * circle;  % inertial frame
 
            plot3(robot(1,:),robot(2,:),robot(3,:),'Color',[mod(zz,2),0,1-mod(zz,2)]) % RGB/identity

            arrow3([0 0 0],[-0.1 0 0])           % X
            arrow3([0 0 0],[0 0.1 0])            % Y  inertial frame
            arrow3([0 0 0],[0 0 0.1])            % Z

            text(-0.1,0,0,'Y') 
            text(0,0.1,0,'X')                    %% corresponding to body frame
            text(0,0,0.1,'Z')

        end

    end
    
    F   = getframe(gcf); 

    writeVideo(video,F);
    
end

%---------functions for solving 'g'------------%

function skew=vector_tilde(vec)

    skew        = zeros(3,3);
    skew(1,2)   = -vec(3);
    skew(1,3)   = vec(2);
    skew(2,1)   = vec(3);
    skew(2,3)   = -vec(1);
    skew(3,1)   = -vec(2);
    skew(3,2)   = vec(1); 

end

function matrix = vector_hat(screw)

    matrix         = zeros(4,4);
    matrix(1:3,1:3)= vector_tilde(screw(1:3)); 
    matrix(1:3,4)  = screw(4:6);

end

function g = piecewise_expmap(x,theta,xi)

    xihat           = vector_hat(xi);
    if theta==0
        g           = eye(4)+x*xihat;

    else
        g           = eye(4)+x*xihat+...
                      ((1-cos(x*theta))/(theta^2))*xihat^2+...
                      ((x*theta-sin(x*theta))/(theta^3))*xihat^3;
    end
end