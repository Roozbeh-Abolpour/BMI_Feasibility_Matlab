function [x,y,g,dt]=InnerLoop(Ls,Bs,bounds,x0,y0,g0,conv_ther,tmax)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
nx=length(x0);ny=length(y0);
nl=length(Ls);nb=length(Bs);
x=x0;y=y0;g=g0;
dt=0;
while dt<=tmax
    tic
    gb=g;xb=x;
    yalmip('clear');
    x=sdpvar(nx,1);g=sdpvar(1,1);    
    F=[x>=lbx,x<=ubx];    
    for i=1:nl
        L=Ls{i}(x,y);
        F=[F,L<=g*eye(length(L))];
    end
    for i=1:nb
        B=Bs{i}(x,y);
        F=[F,B<=g*eye(length(B))];
    end
    op=sdpsettings;op.verbose=0;op.solver='mosek';    
    res=optimize(F,g,op);    
    if res.problem~=0
        x=xb;
        dt=dt+toc;
        return
    end
    g=value(g);x=value(x);        
    if abs(g-gb)<=conv_ther||g<0  
        dt=dt+toc;
        break
    end
    yb=y;
    yalmip('clear');
    y=sdpvar(ny,1);g=sdpvar(1,1);    
    F=[y>=lby,y<=uby];    
    for i=1:nl
        L=Ls{i}(x,y);
        F=[F,L<=g*eye(length(L))];
    end
    for i=1:nb
        B=Bs{i}(x,y);
        F=[F,B<=g*eye(length(B))];
    end
    op=sdpsettings;op.verbose=0;op.solver='mosek';    
    res=optimize(F,g,op);    
    if res.problem~=0
        y=yb;
        dt=dt+toc;
        return
    end
    g=value(g);y=value(y);   
    if abs(g-gb)<=conv_ther||g<0        
        dt=dt+toc;
        break
    end
    dt=dt+toc;
end
end