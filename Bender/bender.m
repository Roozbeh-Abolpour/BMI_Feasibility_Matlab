function [x,y,T]=bender(Ls,Bs,bounds,params)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;
nx=length(lbx);ny=length(lby);
[Lss,Bss]=extendedform(Ls,Bs,nx,ny);
Uss={};Vss={};us={};vs={};
conv_ther=params.conv_ther;
accept_ther=params.accept_ther;
tmax=params.tmax;
x=0.5*(lbx+ubx);
T=0;
while T<=tmax
    tic
    xb=x;
    x=master(Lss,Bss,bounds,Uss,Vss,us,vs);
    if isempty(x)
        x=[];y=[];
        T=T+toc;
        return
    end
    [g,y]=sub(Ls,Bs,bounds,x);
    [Us,Vs,u,v]=dualsub(Lss,Bss,bounds,x);
    if g<=accept_ther||norm(x-xb)<=conv_ther
        T=T+toc;
        return
    end 
    Uss{end+1}=Us;Vss{end+1}=Vs;
    us{end+1}=u;vs{end+1}=v;
    T=T+toc;
end
end