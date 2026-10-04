function [xb,yb,g]=InnerLoop(Ls,Bs,x0,y0)
xb=x0;yb=y0;conv_tol=1e-4;
g=ErrBound(Ls,Bs,xb,yb);
while 1==1    
    gb=g;
    [xb,yb]=LMIPart(Ls,Bs,xb,yb,2);
    [xb,yb]=LMIPart(Ls,Bs,xb,yb,1);   
    g=ErrBound(Ls,Bs,xb,yb);
    if abs(g-gb)<=conv_tol||g<=1e-6
        break
    end
end
end
