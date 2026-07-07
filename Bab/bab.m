function [x,y,T]=bab(Ls,Bs,bounds,params)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
nx=length(lbx);ny=length(lby);
[Lss,Bss]=extendedform(Ls,Bs,nx,ny);
conv_length=params.conv_length;
accept_ther=params.accept_ther;
tmax=params.tmax;
Ns{1}=[lbx ubx];bs=inf;vs=max(ubx-lbx);
x=[];y=[];
T=0;
while ~isempty(Ns)&&T<=tmax
    tic
    indh=find(vs<=conv_length);        
    Ns(indh)=[];bs(indh)=[];vs(indh)=[];
    if isempty(Ns)
        T=T+toc;
        return
    end
    [~,ind]=min(bs);
    Ne=Ns{ind};
    xe=mean(Ne')';    
    mus=zeros(1,nx);    
    for i=1:nx
        mus(i)=min(xe(i)-Ne(i,1),Ne(i,2)-xe(i));
    end
    [~,inde]=max(mus);
    Ne1=Ne;Ne2=Ne;
    Ne1(inde,2)=0.5*(Ne(inde,1)+Ne(inde,2));
    Ne2(inde,1)=0.5*(Ne(inde,1)+Ne(inde,2));
    Ns(ind)=[];bs(ind)=[];vs(ind)=[];
    for i=1:2
        if i==1
            Ne=Ne1;
        else
            Ne=Ne2;
        end
        lbx=Ne(:,1);ubx=Ne(:,2);
        [xe,~,~,gl]=bablb(Lss,Bss,lbx,ubx,lby,uby);
        if isempty(xe)||(~isempty(bs)&&gl>min(bs))
            continue
        end           
        xe1=mean(Ne')';    
        vs(end+1)=max(Ne(:,2)-Ne(:,1));
        Ns{end+1}=Ne;
        [ye1,g1]=babuby(xe1,Ls,Bs,lby,uby);        
        if isempty(g1)
            bs(end+1)=inf;            
            continue
        end 
        if g1<=accept_ther
            x=xe1;y=ye1;
            T=T+toc;
            return
        elseif g1<min(bs)
            x=xe1;y=ye1;
        end
        ye2=ye1;
        [xe2,g2]=babubx(ye2,Ls,Bs,lbx,ubx);
        if isempty(g2)
            bs(end+1)=g1;
            continue
        end                
        bs(end+1)=g2;
        if g2<=accept_ther
            x=xe2;y=ye2;            
            T=T+toc;
            return
        elseif g2<min(bs)
            x=xe2;y=ye2;
        end
    end
    T=T+toc;
end
end