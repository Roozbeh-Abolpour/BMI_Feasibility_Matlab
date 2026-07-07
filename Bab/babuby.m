function [y,g]=babuby(x,Ls,Bs,lby,uby)
ny=length(lby);
yalmip('clear')
y=sdpvar(ny,1);
g=sdpvar(1,1);
F=[y>=lby,y<=uby];
for i=1:length(Ls)
    L=Ls{i}(x,y);
    F=[F,L<=g*eye(length(L))];
end
for i=1:length(Bs)
    B=Bs{i}(x,y);
    F=[F,B<=g*eye(length(B))];
end
ops=sdpsettings;ops.solver='mosek';ops.verbose=0;
sol=optimize(F,g,ops);
if sol.problem==0
    g=value(g);y=value(y);
else
    y=[];g=inf;   
end
end