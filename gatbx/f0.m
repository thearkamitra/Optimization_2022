function [s]=f0(Em1,Em2,L1,L2,L3,Ek0,a)  
 syms Vt p B z y
     zt=[];
     t=[];
     pt=[];
     E=[];
     p1=[];
     Vp1=[];
     Vt1=[];
     t1=[];
     EF = [];
     ts=[];
     e=1.6*10^-19;
     m0=9.1*10^-31;
     c=3*10^8;
     z0=0;
     p0=sqrt(((Ek0*10^3*1.6*10^-19+m0*c^2)^2-m0^2*c^4)/c^2)/(m0*c);
     w=1.8*10^10;
     flag1=0;
     flag2=0;
     flag3=0;
     flag=0;
     pos=0;
     arc=[];
     Eavg=0;
     for j=-360:5:360
         ti=j*pi/(180*w);
        for i=0:1:8000
            t(i+1)=10^-13.*i+j*pi/(180*w);
            if z0<L1
                Vp=e*Em1*10^6/(m0*c)*sin(pi*z0/L1)*sin(w*t(i+1));
            elseif z0<(L1+L2)
                Vp=0;
            elseif z0<(L1+L2+L3)
                Vp=-e*Em2*10^6/(m0*c)*sin(pi*(z0-L1-L2)/L3)*sin(w*t(i+1)+a);
            end
            y=sqrt(1+p0^2);
            if p0>=0
                B=sqrt(p0^2/(1+p0^2));
            else
                B=-sqrt(p0^2/(1+p0^2));
            end
            Vt=B*c;
            z=z0+Vt*10^-13;
            z0=z;
            zt(i+1)=z0;
            pt(i+1)=p0;
            if z0<0
                flag1=flag1+1;
%                 figure(1)
%                 plot(t,zt,'-r')
%                 hold on;
                break
            elseif z0>(L1+L2+L3)
                flag2=flag2+1;
                E(flag2)=(y-1)*m0*c^2;
                Eavg=mean(E);
                p1(flag2)=p0*(m0*c);
                Vp1(flag2)=Vp;
                Vt1(flag2)=Vt;
                t1(flag2)=t(i+1);
                arc(flag2)=j;
%                 figure(1)
%                 plot(t,zt,'-b')
%                 hold on;
                if((t(i+1)-ti)>10^-9)
                    flag3=flag3+1;
                end
                break
            else 
                p=p0+Vp*10^-13;
                p0=p;
            end
        end
    flag=-flag1+flag2; 
    zt=[];
    pt=[];
    t=[];
    z0=0;
    Vp=0;
    Vt=0;
    z=0;
    p=0;
    p0=sqrt(((Ek0*10^3*1.6*10^-19+m0*c^2)^2-m0^2*c^4)/c^2)/(m0*c);
        end
     if(flag3==0&& Eavg>=1*10^6*1.6*10^-19)
         s=flag2/(flag1+flag2);
         s(isnan(s)==1) = 0;
     else
         s=0;
     end

end
    