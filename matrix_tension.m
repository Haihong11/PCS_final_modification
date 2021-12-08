function T = matrix_tension(t,cavo)
global gv

t_tens      =gv.t_tens;
% T11        =gv.T11;
T13         =gv.T13;
T21         =gv.T21;
% T31        =gv.T31;
T32         =gv.T32;
% T33        =gv.T33;

if cavo==11
%     T=interp1(t_tens,T11,t);
    T=0;
end

if cavo==12
    T=0;
end

if cavo==13
    T=interp1(t_tens,T13,t);
%     T=0;
end