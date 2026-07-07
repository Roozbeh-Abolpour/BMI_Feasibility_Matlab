function [x,y,T]=bpath(Ls,Bs,bounds,params)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
conv_ther=params.conv_ther;
path_radius=params.path_radius;
init_violation=params.init_violation;
accept_ther=params.accept_ther;
tmax=params.tmax;
T=0;
tic
nx=length(lbx);ny=length(lby);
y=0.5*(lby+uby);gb=init_violation;
yalmip('clear');
x=sdpvar(nx,1);
g=sdpvar(1,1);
F=[x>=lbx,x<=ubx];
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
    g=value(g);x=value(x);
else
    x=[];y=[];
    return
end
T=T+toc;
while T<=tmax
    tic
    if g<=accept_ther
        T=T+toc;
        return
    end
    if abs(g-gb)<=conv_ther
        T=T+toc;
        return
    end
    gb=g;   
    yalmip('clear');
    dx=sdpvar(nx,1);
    dy=sdpvar(ny,1);
    g=sdpvar(1,1);
    F=[x+dx>=lbx,x+dx<=ubx];
    F=[F,y+dy>=lby,y+dy<=uby];
    for i=1:length(Ls)
        L=Ls{i}(x+dx,y+dy);
        F=[F,L<=g*eye(length(L))];
    end
    for i=1:length(Bs)
        B=Bs{i}(x,y)+Bs{i}(dx,y)+Bs{i}(x,dy);
        F=[F,B<=g*eye(length(B))];
    end
    F=[F,[-eye(nx) dx;dx' -path_radius]<=0,...
        [-eye(ny) dy;dy' -path_radius]<=0];
    ops=sdpsettings();ops.solver='mosek';ops.verbose=0;
    sol=optimize(F,g,ops);
    if sol.problem==0
        x=x+value(dx);
        y=y+value(dy);
        g=value(g);        
    else
        T=T+toc;
        x=[];y=[];
        return
    end    
    T=T+toc;
end
end