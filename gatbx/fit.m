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
     Ek0=14.494*10^3*1.6*10^-19;
     a0=0;
     e=1.6*10^-19;
     m0=9.1*10^-31;
     Em=25.491*10^6;
     L=0.0276;
     c=3*10^8;
    z0=0;
     p0=sqrt(((Ek0+m0*c^2)^2-m0^2*c^4)/c^2)/(m0*c);
     w=1.8*10^10;
     flag1=0;
     flag2=0;
     flag=0;
     count = 0;
     sum = 0;
     qualified = [];
     for j=-360:5:360
        for i=0:1:5000
            t(i+1)=10^-13.*i+j*pi/(180*w);
            Vp=e*Em/(m0*c)*sin(pi*z0/L)*sin(w*t(i+1));
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
                break
            elseif z0>L
                flag2=flag2+1;
                E(flag2)=(y-1)*m0*c^2/e;
                p1(flag2)=p0*(m0*c);
                Vp1(flag2)=Vp;
                Vt1(flag2)=Vt;
                t1(flag2)=t(i+1);
                break
            else 
                p=p0+Vp*10^-13;
                p0=p;
    end
  end
    flag=-flag1+flag2; 
    figure(1)
    plot(t,zt,'-b')
    hold on;
%     figure(2)
%     plot(t,pt,'-g');
%     hold on;
    zt=[];
    pt=[];
    t=[];
    z0=0;
    Vp=0;
    Vt=0;
    z=0;
    p=0;
    p0=sqrt(((Ek0+m0*c^2)^2-m0^2*c^4)/c^2)/(m0*c);
    end
     flag1;
     flag2;
     flag;
     E;
     p1;
     Vp1;
     Vt1;
     s=flag2/(flag1+flag2);
      tmin=min(t1(:));
       gra = gradient(t1);
       che = gradient(gra);
     for k=1:1:flag2
         if t1(k)<tmin+2*10^-11
             qualified(count+1) = k;
             count = count+1;
         end
     end;
     for k = 1:1:count
         EF(k) = E(k);
     end
       dE=std(EF);
       gra = gradient(t1);
       che = gradient(gra);
       edge = 3*mean(abs(gra));
       pos=int8(0)
     for i=1:length(t1)
         if che(i)<0
             pos=i+1;
             break
         end
     end
     for i=(pos+1):1:length(t1)
         if t1(i)<=t1(pos)
             ts(i-pos+1)=t1(i);
         else
             break;
         end
     end
     dT=std(ts);
     t1s = [];
     gras = [];
     ches = [];
     for i = 1:1:length(t1)
         t1s(i,1) = i;
         t1s(i,2) = t1(i);
     end
     for i = 1:1:length(gra)
         gras(i,1) = i;
         gras(i,2) = gra(i);
     end
     for i = 1:1:length(che)
         ches(i,1) = i;
         ches(i,2) = che(i);
     end
     figure(3)
     plot(t1s(:,1),t1s(:,2),'r*')
     title('t1')
     figure(4)
     plot(gras(:,1),gras(:,2),'r*')
     title('gra')
     figure(5)
     plot(ches(:,1),ches(:,2),'r*')
     title('che')
