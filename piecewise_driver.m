clear
clc
format long

global gv
warning off

tic

disp('Pre-processing')

%---------Get geometrical parameters----------%
 
E         = 110e3;                          % [Pa] Young modulus
mu        = 300;                            % [Pa*s] Shear viscosity modulus
Poi       = 0.45;                           % [-]  Poisson ratio
G         = E/(2*(1+Poi));                  % [Pa] Shearing modulus
L         = [0.09,0.07,0.04];           % [m]  Array consist of length of soft robot of each section
num_piece = 3;                              % Number of piece 
num_disc  = [10,8,5];                        % Array consist of number of disc for each section.
%---------coordinate of each disc-----------%
X1        = linspace(0,L(1),num_disc(1));   
X2        = X1(end) + linspace(0,L(2),num_disc(2));
X3        = X2(end) + linspace(0,L(3),num_disc(3));
X1(end)   = [ ];
X2(end)   = [ ];
X         = [X1,X2,X3];
dX        = L./(num_disc-1);                % The spacing between the discs
Rmax      = 0.01;                         % [m] Radius of base
Rmin      = 0.005;                         % [m] Radius of end
ro_arm    = 2000;                           % [kg/m^3]
Gra       = 0*[0;0;0;0;0;-9.81];              % [m/s^2] gravity acceleration twist in the spatial frame

xi_star   = [0;0;0;1;0;0];                  % the reference configuration 

%--------------Cable tension------------%

T           = [7,0,0,0];                    % [N]  Tension of the four cables (Scalar)
% d           = 7e-3;
d           = [1e-2,9e-3,7e-3];             % Distance between the cable and midline from the first section to last one.

%--------------Numerical setting--------------%

time        = 1;                            % [s]
nsol        = time*10^2+1;                  % Number of solution 
tspan       = linspace(0,time,nsol);        % [0 time];    

gv.ro_arm      = ro_arm;
gv.E           = E;
gv.mu          = mu;
gv.G           = G;
gv.Gra         = Gra;
gv.L           = L;
gv.X           = X;
gv.xi_star     = xi_star;
gv.nsol        = nsol;
gv.num_disc    = num_disc;
gv.num_piece   = num_piece;
gv.dX          = dX;
gv.time        = time;
gv.tspan       = tspan;
gv.T           = T;
gv.d           = d;
gv.Rmax        = Rmax;
gv.Rmin        = Rmin;

disp('Time-advancing')

options        = odeset('RelTol',1e-4);

%------------------Initial conditions---------------%

xi_0           = [0;0;0;1;0;0];  % given variable of the joint, column vector
xidot_0        = [0;0;0;0;0;0];
    
ini_cond       = [repmat(xi_0',[1,num_piece]) repmat(xidot_0',[1,num_piece])];  % row array

%------------------RK-45------------------%

[t,z]          = ode45(@piecewise_CBA,tspan,ini_cond,options);

toc

disp('Post-processing')

%--------------Define 6 strain components for each section----------%

% s         = zeros(nsol,num_piece);
% v         = zeros(nsol,num_piece);
% k         = zeros(nsol,num_piece);
% q         = zeros(nsol,num_piece);
% p         = zeros(nsol,num_piece);
% r         = zeros(nsol,num_piece);
% 
% for ii=1:num_piece 
% 
%     s(:,ii)   = z(:,6*(ii-1)+1); 
%     v(:,ii)   = z(:,6*(ii-1)+2);
%     k(:,ii)   = z(:,6*(ii-1)+3);
%     q(:,ii)   = z(:,6*(ii-1)+4);
%     p(:,ii)   = z(:,6*(ii-1)+5);
%     r(:,ii)   = z(:,6*(ii-1)+6);  
%     
% end
%  
% %--------------Plot--------------%
% 
% for ii=1:num_piece
%  
%     figure
%     plot(t,s(:,ii))
%     grid on
%     title(strcat('torsion of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('s [1/m]')
%     print('-djpeg')
% 
%     figure
%     plot(t,v(:,ii))
%     grid on
%     title(strcat('curvature on y of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('v [1/m]')
%     print('-djpeg')
% 
%     figure
%     plot(t,k(:,ii))
%     grid on
%     title(strcat('curvature on z of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('k [1/m]')
%     print('-djpeg')
% 
%     figure
%     plot(t,q(:,ii))
%     grid on
%     title(strcat('longitudinal strain of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('q [1]')
%     print('-djpeg')
% 
%     figure
%     plot(t,p(:,ii))
%     grid on
%     title(strcat('tras y strain of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('p [1]')
%     print('-djpeg')
% 
%     figure
%     plot(t,r(:,ii))
%     grid on
%     title(strcat('tras z strain of piece',num2str(ii)))
%     xlabel('t [s]')
%     ylabel('r [1]')
%     print('-djpeg')
%    
% end