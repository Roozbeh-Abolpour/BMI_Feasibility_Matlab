function [x,y,T]=mine(Ls,Bs,bounds,params)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
conv_ther=params.conv_ther;
accept_ther=params.accept_ther;
sep_cap=params.sep_cap;
tmax=params.tmax;
T=0;
tic
x=0.5*(lbx+ubx);y=0.5*(lby+uby);g=0;
for i=1:length(Ls)
    L=Ls{i}(x,y);
    mu=max(real(eig(L)));
    if mu>0
        g=max(g,mu);
    end
end
for i=1:length(Bs)
    B=Bs{i}(x,y);
    mu=max(real(eig(B)));
    if mu>0
        g=max(g,mu);
    end
end
T=T+toc;
while T<=tmax    
    gb=g;
    [xb,yb,g,dt]=InnerLoop(Ls,Bs,bounds,x,y,g,conv_ther,tmax-T);        
    T=T+dt;
    x=xb;y=yb;
    if g<=accept_ther||abs(g-gb)<=conv_ther        
        break
    end
    tic
    [Ps,Qs]=SepPlane(Ls,Bs,xb,yb,sep_cap);     
    T=T+toc;
    if isempty(Ps)
        break
    end
    str='Ls{end+1}=@(x,y) ';
    for i=1:length(Ls)
        str=[str sprintf('trace(Ps{%d}*Ls{%d}(x,y))',i,i) '+'];
    end
    for i=1:length(Bs)
        str=[str sprintf('trace(Qs{%d}*Bs{%d}(x,y))',i,i) '+'];
    end
    str=[str(1:end-1) ';'];
    eval(str);             
end
end