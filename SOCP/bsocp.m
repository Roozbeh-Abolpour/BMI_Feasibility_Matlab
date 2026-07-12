function [x,y,T]=bsocp(Ls,Bs,bounds,params)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
nx=length(lbx);ny=length(lby);
T=0;
tic
[~,Bss]=extendedform(Ls,Bs,nx,ny);
yalmip('clear');
x=sdpvar(nx,1);y=sdpvar(ny,1);
Z=sdpvar(nx+ny,nx+ny,'symmetric');
g=sdpvar(1,1);
F=[x>=lbx,x<=ubx,y>=lby,y<=uby];
for i=1:length(Ls)
    L=Ls{i}(x,y);
    F=[F,L<=g*eye(length(L))];
end
for i=1:length(Bss)
    B=Bss{i}{1};
    for j=1:nx
        B=B+Bss{i}{1+j}*x(j);
    end
    for j=1:ny
        B=B+Bss{i}{1+nx+j}*y(j);
    end    
    for j=1:nx
        for k=1:ny            
            B=B+Bss{i}{1+nx+ny+(j-1)*ny+k}*Z(j,nx+k);
        end
    end
    F=[F,B<=g*eye(length(B))];
end
z=[x;y];
M=[-1 z';z -Z];
for i=1:1+nx+ny
    F=[F,M(i,i)<=0];
    for j=1:i-1
        a=M(i,i);b=M(i,j);c=M(j,j);
        F=[F,c<=0];
        F=[F,cone([2*b; a-c],-a-c)];
    end
end
ops=sdpsettings;ops.solver='mosek';ops.verbose=0;
sol=optimize(F,g,ops);
if sol.problem==0
    x=value(x);y=value(y);
else
    x=[];y=[];    
end
T=T+toc;
end
