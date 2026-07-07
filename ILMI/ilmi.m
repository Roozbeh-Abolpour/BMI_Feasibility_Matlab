function [x,y,T]=ilmi(Ls,Bs,bounds,params)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
conv_ther=params.conv_ther;
accept_ther=params.accept_ther;
tmax=params.tmax;
init_violation=params.init_violation;
nx=length(lbx);ny=length(lby);
y=0.5*(lby+uby);gb=init_violation;
T=0;
while T<=tmax
    tic
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
        T=T+toc;
        return
    end
    yalmip('clear');
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
        x=[];y=[];
        T=T+toc;
        return
    end
    if g<=accept_ther
        T=T+toc;
        return
    end
    if abs(g-gb)<=conv_ther
        T=T+toc;
        return
    end
    gb=g;
    T=T+toc;
end
end