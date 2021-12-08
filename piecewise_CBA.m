function dz = piecewise_CBA(t,z)

t

global gv

L           = gv.L;
E           = gv.E;
mu          = gv.mu;
G           = gv.G;
ro_arm      = gv.ro_arm;
Gra         = gv.Gra;
dX          = gv.dX;      %    dX =   L/(num_disc-1);
X           = gv.X;       %     X =   linspace(0,L,num_disc)
T           = gv.T;       %    Array consist of cable tension.
d           = gv.d;
num_piece   = gv.num_piece;
num_disc    = gv.num_disc;
Rmax        = gv.Rmax;
Rmin        = gv.Rmin;

% Code for getting parameters. 

xi_star     = [0;0;0;1;0;0];             %--Compute the internal elastic force

g_r         = [1  0 0 0; 
               0  1 0 0; 
               0  0 1 0;
               0  0 0 1];                %  the body frame relative to the spatial frame

eta_pre     = zeros(6,1);                % velocity of the first disc
g_pre       = eye(4); 

GIM         = zeros(6*num_piece);        % Generalized inertia matrix
GCM1        = zeros(6*num_piece);        % Generalized Coriolis matrix 1
GCM2        = zeros(6*num_piece);        % Generalized Coriolis matrix 2
GM          = zeros(6*num_piece,6);      % Gravitational matrix

% Current system state

Xi          = z(1:6*num_piece);                       % strain twist 
Xidot       = z(6*num_piece+1:12*num_piece);          % derivative of strain

%------------Get the configuration and velocity of all discs-----------%

for i = 1:num_piece

    xin          = Xi(6*i-5:6*i);                 % column vector row array
    xidotn       = Xidot(6*i-5:6*i);
    kn           = xin(1:3);
    thetan       = sqrt(kn' * kn);

    if i == 1

        disc_index = 1;

    else

        disc_index = sum(num_disc(1:i-1)) + 1;

    end
      
    num_disc_i   = num_disc(i);                
    
    for ii  = 1:num_disc_i
   
        invAdjg      = piecewise_invAdjoint(X(ii),thetan,xin);
  
        intdAdjg     = piecewise_ADJ(X(ii),thetan,xin);

        eta(:,(disc_index-1)+ii) = invAdjg * (eta_pre+intdAdjg * xidotn);

        g(:,4*(disc_index-1)+4 * ii-3:4*(disc_index-1) + 4 * ii)...
                         = g_pre * piecewise_expmap(X(ii),thetan,xin); 
        
    end
    
    invAdjg_next      = piecewise_invAdjoint(X(num_disc_i),thetan,xin);
    
    intdAdjg_next     = piecewise_ADJ(X(num_disc_i),thetan,xin);
    
    eta_pre           = invAdjg_next * (eta_pre+intdAdjg_next * xidotn);
    
    g_pre             = g_pre * piecewise_expmap(X(num_disc_i),thetan,xin);

end

%-----Compute Jacobian before actual using it (body Jacobian)---------%.

J_i_list = [];                     % Save Jacobian matrix of all discs

J_pre    = zeros(6,6*num_piece);   % initial Jacobian matrix

for i = 1:num_piece 
    
    xin          = Xi(6*i-5:6*i); 
    kn           = xin(1:3);
    thetan       = sqrt(kn' * kn);
    
    num_disc_i   = num_disc(i);
    
    for j = 1:num_disc_i       
        
        invAdjg        = piecewise_invAdjoint(X(j),thetan,xin);
        
        intdAdjg       = piecewise_ADJ(X(j),thetan,xin);      %  T_g(X).
        
        J_i =  J_pre;                           
        
        J_i(:, 6*i-5:6*i) = intdAdjg; 
        
        J_i = invAdjg * J_i;
        
        if j == num_disc_i                

            J_pre = J_i;                 

        end

        J_i_list = [J_i_list; J_i];                         %  6*piece*disc x 6*piece.

    end

end

% Compute the derivative of Jacobian before actual using it.  ( Note the Integral item)

dotJ_i_list = [ ];                    % Save derivative of Jacobian matrix of all discs

dotJ_pre    = zeros(6,6*num_piece); 

for i = 1:num_piece
    
    xin          = Xi(6*i-5:6*i); 
    kn           = xin(1:3);
    thetan       = sqrt(kn' * kn);
          
    if i == 1

        disc_index = 1;

    else

        disc_index = sum(num_disc(1:i-1)) + 1;

    end
    
    %----Analyse 'eta_pre'----%
    
    ini_disc     = 2;

    ii           = ini_disc-1;

    if i ~= 1

        eta_pre = eta(:,(disc_index-1)+ii);

    else

        eta_pre = zeros(6,1);

    end      
    
    num_disc_i   = num_disc(i);
    
    for ii = 1:num_disc_i
    
      dotJ_i         = dotJ_pre;

      invAdjg        = piecewise_invAdjoint(X(ii),thetan,xin);  %  Ad^{-1}_g(X)
 
      intdAdjg       = piecewise_ADJ(X(ii),thetan,xin);

      dotJ_i(:, 6*i-5:6*i) = matrix_adj(eta_pre) * intdAdjg;        

      dotJ_i               = invAdjg * dotJ_i;

      if  ii == num_disc_i               

        dotJ_pre = dotJ_i;                 

      end

      dotJ_i_list          = [dotJ_i_list; dotJ_i];     %  6 * piece * disc x 6 * piece.

    end

end

%--------Compute the variable M,Eps,Ips, e.g.,the variable mass matrix of discs etc.------%

Eps_list = [ ];

Ips_list = [ ];

M_list   = [ ];

R_d = (Rmax-Rmin) / sum(L);                     % Radius reduction of the circular section for unit length X.

for i = 1:num_piece
    
    num_disc_i   = num_disc(i);
    
    for ii = 1:num_disc_i
        
        R_ii       = Rmax - R_d * (ii-1) * dX(i);   
        
        A          = pi*R_ii^2;                         % [m^2] Cross-section area  
        J          = pi*R_ii^4/4;                       % [m^4] The second moment of the area, e.g.,'J_y' or 'J_z' 
        I          = pi*R_ii^4/2;                       % [m^4] The polar meomet of the area e.g., I_x = J_y + J_z

        Eps        = diag([G*I E*J E*J E*A G*A G*A]);  % stifness matrix
        
        Ips        = mu * diag([I 3*J 3*J 3*A A A]);   % viscosity matrix

        M          = ro_arm * diag([I J J A A A]);     % screw inertia matrix of the disc.

        if  ii == num_disc_i

            Rmax   = R_ii;

        end

        Eps_list = [Eps_list;Eps];

        Ips_list = [Ips_list;Ips];

        M_list   = [M_list;M];

    end
    
end

%-----------------Composite Body Algorithm (CBA)--------------------------%

for m = 1:num_piece 

    for n = 1:num_piece 
    
        M_mn   = zeros(6);
        C1_mn  = zeros(6);
        C2_mn  = zeros(6);
        
        for i = max(m,n):num_piece
            
            xin          = Xi(6*i-5:6*i); 
            xidotn       = Xidot(6*i-5:6*i);
            kn           = xin(1:3);
            thetan       = sqrt(kn' * kn);
            M_mn_i       = zeros(6);
            C1_mn_i      = zeros(6);
            C2_mn_i      = zeros(6);
            
            %-------Calculate the first disc of certain section-------%
                
            ini_disc   = 2;

            ii         = ini_disc-1;

            if i == 1
                
                disc_index = ii;
                
            else

                disc_index = sum(num_disc(1:i-1)) + ii;
  
            end

            J_i        = J_i_list(6*disc_index-5:6*disc_index,:);

            dotJ_i     = dotJ_i_list(6*disc_index-5:6*disc_index,:);

            Sm         = J_i(:,6*m-5:6*m);

            Sn         = J_i(:,6*n-5:6*n);

            dotSn      = dotJ_i(:,6*n-5:6*n);
            
            %--------------Variable mass of disc----------%
            
            M_ii       = M_list(6*disc_index-5:6*disc_index,:);
            
            %--------------------Mass--------------------%

            M_p        = Sm' * M_ii * Sn;

            %--------------------GCM1--------------------%
            
            if i ~= 1
                
                eta_pre = eta(:,(disc_index-1)+ii);

            else
 
                eta_pre = zeros(6,1);

            end

            C_p1       = Sm' * matrix_coadj(eta_pre) * M_ii * Sn;

            %--------------------GCM2-------------------%

            C_p2       = Sm'* M_ii * dotSn;

            %-------------------------------------------%

            num_disc_i = num_disc(i);

            for ii = ini_disc:num_disc_i

                invAdjg    = piecewise_invAdjoint(X(ii),thetan,xin);

                intdAdjg   = piecewise_ADJ(X(ii),thetan,xin);

                eta_disc   = invAdjg*(eta_pre+intdAdjg * xidotn);

                if i == 1
                    
                    disc_index = ii;
                    
                else
                    
                    disc_index = sum(num_disc(1:i-1)) + ii;
                    
                end

                J_i        = J_i_list(6*disc_index-5:6*disc_index,:);

                dotJ_i     = dotJ_i_list(6*disc_index-5:6*disc_index,:);

                Sm         = J_i(:,6*m-5:6*m); 

                Sn         = J_i(:,6*n-5:6*n);

                dotSn      = dotJ_i(:,6*n-5:6*n);
                
                M_ii       = M_list(6*disc_index-5:6*disc_index,:);

                M_mn_i     = M_mn_i + (M_p + Sm' * M_ii * Sn) /2 * dX(i);                              % equation 31 

                C1_mn_i    = C1_mn_i + (C_p1 + Sm' * matrix_coadj(eta_disc) * M_ii * Sn) /2 * dX(i);   % equation 32

                C2_mn_i    = C2_mn_i + (C_p2 + Sm'* M_ii * dotSn) /2 * dX(i);                          % equation 33
 
                M_p        = Sm' * M_ii * Sn;

                C_p1       = Sm' * matrix_coadj(eta_disc) * M_ii * Sn;

                C_p2       = Sm' * M_ii * dotSn;
                
             end

           %------Sum again after integrating-----%
            
            M_mn  = M_mn  + M_mn_i;
            C1_mn = C1_mn + C1_mn_i;
            C2_mn = C2_mn + C2_mn_i;

        end

        %  Fill in the integrated block.
        
        GIM(6*m-5:6*m, 6*n-5:6*n)  = M_mn;   
        GCM1(6*m-5:6*m, 6*n-5:6*n) = C1_mn;
        GCM2(6*m-5:6*m, 6*n-5:6*n) = C2_mn;

    end       

end

%-------------------------Compute gravitational matrix-------------------------% 

for m = 1:num_piece         % 'm' represents the number of the section

    G_m         = zeros(6);
     
        for i = m:num_piece  

            xin          = Xi(6*i-5:6*i);                      
            kn           = xin(1:3); 
            thetan       = sqrt(kn'*kn);
            G_m_i        = zeros(6);
           
            %-------Calculate the first disc-------%

            ini_disc   = 2;

            ii         = ini_disc-1;

            if i == 1 
                
                disc_index = ii;
                
            else
                
                disc_index = sum(num_disc(1:i-1)) + ii;
                
            end

            J_i          = J_i_list(6*disc_index-5:6*disc_index,:);

            Sm           = J_i(:,6*m-5:6*m);
            
            M_ii         = M_list(6*disc_index-5:6*disc_index,:);
            
            %--------Analyse 'g_pre'--------%
            
            if i ~= 1
          
                g_pre      = g(:,4*(disc_index-1)+4*ii-3:4*(disc_index-1) + 4*ii);  
 
            else
             
                g_pre      = eye(4);
                
            end
            
            G_p          = Sm' * M_ii * matrix_Adjoint(inv(g_pre));
                                              
            num_disc_i   = num_disc(i);

            for ii = ini_disc:num_disc_i

                g_disc     = g_pre * piecewise_expmap(X(ii),thetan,xin);

                if i == 1 
                
                    disc_index = ii;
                
                else
                
                    disc_index = sum(num_disc(1:i-1)) + ii;
                
                end

                J_i        = J_i_list(6*disc_index-5:6*disc_index,:);

                Sm         = J_i(:,6*m-5:6*m);
                
                M_ii       = M_list(6*disc_index-5:6*disc_index,:);
    
                G_m_i      = G_m_i + (G_p + Sm' * M_ii * matrix_Adjoint(inv(g_disc)))/2 * dX(i); 

                G_p        = Sm' * M_ii * matrix_Adjoint(inv(g_disc));

            end

            G_m   =  G_m + G_m_i;

        end
 
        GM(6*m-5:6*m,:)    = G_m;

end

%------------------Internal Load (Fi_n)------------------------% 

Fin_list = [ ];

for i = 1:num_piece
    
    xin         = Xi(6*i-5:6*i,:); 

    xidotn      = Xidot(6*i-5:6*i,:); 
   
    Eps_s       = zeros(6);  

    Ips_s       = zeros(6);
    
    num_disc_i  = num_disc(i);
    
    for ii = 1:num_disc_i
        
        if i == 1 
                
            disc_index = ii;
            
        else
              
            disc_index = sum(num_disc(1:i-1)) + ii;
             
        end
        
        Eps_ii         = Eps_list(6*disc_index-5:6*disc_index,:);
        
        Ips_ii         = Ips_list(6*disc_index-5:6*disc_index,:);
        
        Eps_s          = Eps_s + Eps_ii;             
        
        Ips_s          = Ips_s + Ips_ii;       
        
    end
    
    Eps_i       = Eps_s / num_disc_i;         %  Average stifness matrix
    
    Ips_i       = Ips_s / num_disc_i;         %  Average viscosity matrix
    
    Fin         = Eps_i * (xin-xi_star)+Ips_i * xidotn;

    Fi_n        = -L(i) * Fin;

    Fin_list    = [Fin_list;Fi_n]; 

end


%-----------------------Each section driven by cables------------------------%

H         = zeros(6*num_piece,4);          % Under-actuated

Q_pre     = zeros(6,4);

for i = num_piece:-1:1 

    pi1        = [0;d(i);0];           % [m]  cable midline distance i1 (Positive direction of y axis)

    pi2        = [0;0;d(i)];           % [m]  cable midline distance i2 (Positive direction of z axis)

    pi3        = [0;-d(i);0];          % [m]  cable midline distance i3 (Negative direction of y axis)

    pi4        = [0;0;-d(i)];          % [m]  cable midline distance i4 (Negative direction of z axis)

    g_ci1 = [eye(3) pi1;
            zeros(1,3) 1];

    g_ci2 = [eye(3) pi2;
            zeros(1,3) 1];

    g_ci3 = [eye(3) pi3;
            zeros(1,3) 1];

    g_ci4 =  [eye(3) pi4;
            zeros(1,3) 1];
     
    Q_i   = [matrix_coAdjoint(g_ci1) * [0 0 0 -1 0 0]',matrix_coAdjoint(g_ci2) * [0 0 0 -1 0 0]',...
             matrix_coAdjoint(g_ci3) * [0 0 0 -1 0 0]',matrix_coAdjoint(g_ci4) * [0 0 0 -1 0 0]'];   % 6 x 4
         
    Q     = Q_pre + Q_i;
 
    H(6*i-5:6*i,:) = L(i) * Q;
 
    Q_pre = Q;
 
end

%--------------State equation-----------------%

dotz1         = Xidot;
dotz2         = GIM^-1 * (H*T'+Fin_list+GM*matrix_Adjoint(g_r^-1)*Gra-(GCM1+GCM2)*Xidot);
dz            = [dotz1;dotz2];

Tau           = H*T';

gv.Tau        = Tau;
gv.GIM        = GIM;
gv.GM         = GM;
gv.Fin_list   = Fin_list;
gv.GCM1       = GCM1;
gv.GCM2       = GCM2;

end