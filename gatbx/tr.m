function Chrom=tr(Xcluster)
[~,wow]=size(Xcluster); ObjV1t1=0;ObjV1t2=0;ObjV1t3=0;ObjV1t4=0;
a1=[];a2=[];a3=[];a4=[];b1=[];b2=[];b3=[];b4=[];c1=[];c2=[];c3=[];c4=[];d1=[];d2=[];d3=[];d4=[];e1=[];e2=[];e3=[];e4=[];f1=[];f2=[];f3=[];f4=[];g1=[];g2=[];g3=[];g4=[];
if wow>=1
    [o,~]=size(Xcluster{1,1});u=Xcluster{1,1};
end
if wow>=2
    [p,~]=size(Xcluster{1,2});v=Xcluster{1,2};
end
if wow>=3
    [q,~]=size(Xcluster{1,3});w=Xcluster{1,3};
end
if wow>=4
    [h,~]=size(Xcluster{1,4});k=Xcluster{1,4};
end
if wow>=1
    for j=1:o
        tempa=tentotwo1(u(j,1));
        a1(j,:)=tempa;
    end
    for j=1:o
        tempa=tentotwo1(u(j,2));
        b1(j,:)=tempa;
    end
    for j=1:o
        tempa=tentotwo2(u(j,3));
        c1(j,:)=tempa;
    end
    for j=1:o
        tempa=tentotwo2(u(j,4));
        d1(j,:)=tempa;
    end
    for j=1:o
        tempa=tentotwo2(u(j,5));
        e1(j,:)=tempa;
    end
    for j=1:o
        tempa=tentotwo3(u(j,6));
        f1(j,:)=tempa;
    end
    for j=1:o
        if u(j,7)>1
        tempa=tentotwo4(u(j,7));
        else
        tempa=tentotwo2(u(j,7));
        end
        g1(j,:)=tempa;
    end
end

if wow>=2
    for j=1:p
        tempa=tentotwo1(v(j,1));
        a2(j,:)=tempa;
    end
    for j=1:p
        tempa=tentotwo1(v(j,2));
        b2(j,:)=tempa;
    end
    for j=1:p
        tempa=tentotwo2(v(j,3));
        c2(j,:)=tempa;
    end
    for j=1:p
        tempa=tentotwo2(v(j,4));
        d2(j,:)=tempa;
    end
    for j=1:p
        tempa=tentotwo2(v(j,5));
        e2(j,:)=tempa;
    end
    for j=1:p
        tempa=tentotwo3(v(j,6));
        f2(j,:)=tempa;
    end
    for j=1:p
        if v(j,7)>1
        tempa=tentotwo4(v(j,7));
       
        else
        tempa=tentotwo2(v(j,7));
        end
        g2(j,:)=tempa;
    end
        
end
if wow>=3
    for j=1:q
        tempa=tentotwo1(w(j,1));
        a3(j,:)=tempa;
    end
    for j=1:q
        tempa=tentotwo1(w(j,2));
        b3(j,:)=tempa;
    end
    for j=1:q
        tempa=tentotwo2(w(j,3));
        c3(j,:)=tempa;
    end
    for j=1:q
        tempa=tentotwo2(w(j,4));
        d3(j,:)=tempa;
    end
    for j=1:q
        tempa=tentotwo2(w(j,5));
        e3(j,:)=tempa;
    end
    for j=1:q
        tempa=tentotwo3(w(j,6));
        f3(j,:)=tempa;
    end
    for j=1:q
        if w(j,7)>1
        tempa=tentotwo4(w(j,7));
        else
        tempa=tentotwo2(w(j,7));
        end
        g3(j,:)=tempa;
    end
        
end
if wow>=4
    for j=1:h
        tempa=tentotwo1(k(j,1));
        a4(j,:)=tempa;
    end
    for j=1:h
        tempa=tentotwo1(k(j,2));
        b4(j,:)=tempa;
    end
    for j=1:h
        tempa=tentotwo2(k(j,3));
        c4(j,:)=tempa;
    end
    for j=1:h
        tempa=tentotwo2(k(j,4));
        d4(j,:)=tempa;
    end
    for j=1:h
        tempa=tentotwo2(k(j,5));
        e4(j,:)=tempa;
    end
    for j=1:h
        tempa=tentotwo3(k(j,6));
        f4(j,:)=tempa;
    end
    for j=1:h
        if k(j,7)>1
        tempa=tentotwo4(k(j,7));
               else
        tempa=tentotwo2(k(j,7));
         end
        g4(j,:)=tempa;
    end
        
end
Chrom=[a1,b1,c1,d1,e1,f1,g1;a2,b2,c2,d2,e2,f2,g2;a3,b3,c3,d3,e3,f3,g3;a4,b4,c4,d4,e4,f4,g4];
end


function x=tentotwo1(n)
m=fix(n);
t=n-m;
x=[];b=[];a=[];
if m>0
    for i=1:10                      % 处理整数
        a(11-i)=rem(m,2);
        m=fix(m/2);
    end
end
if  t>=0        %处理小数点后的位数，乘2取整法 ，当乘2变为整数后结束
    for i=1:5
        b(i)=fix(t*2);
        t=2*t-fix(2*t);
    end
end
x=[a,b];
end
function x=tentotwo2(n)
m=fix(n);
t=n-m;
x=[];b=[];a=[];
if m>0
    for i=1:1                      % 处理整数
        a(2-i)=rem(m,2);
        m=fix(m/2);
    end
end
if  t>=0        %处理小数点后的位数，乘2取整法 ，当乘2变为整数后结束
    for i=1:10
        b(i)=fix(t*2);
        t=2*t-fix(2*t);
    end
end
x=[a,b];
end
function x=tentotwo3(n)
m=fix(n);
t=n-m;
x=[];b=[];a=[];
if m>0
    for i=1:6                      % 处理整数
        a(7-i)=rem(m,2);
        m=fix(m/2);
    end
end
if  t>=0        %处理小数点后的位数，乘2取整法 ，当乘2变为整数后结束
    for i=1:4
        b(i)=fix(t*2);
        t=2*t-fix(2*t);
    end
end
x=[a,b];
end
% function x=tentotwo4(n)
% m=fix(n);
% t=n-m;
% x=[];b=[];a=[];
% if m>0
%     for i=1:10                      % 处理整数
%         a(11-i)=rem(m,2);
%         m=fix(m/2);
%     end
% end
% if  t>=0        %处理小数点后的位数，乘2取整法 ，当乘2变为整数后结束
%     for i=1:5
%         b(i)=fix(t*2);
%         t=2*t-fix(2*t);
%     end
% end
% x=[a,b];
% end
function x=tentotwo4(n)
m=fix(n);
t=n-m;
x=[];b=[];a=[];
if m>0
    for i=1:2                      % 处理整数
        a(3-i)=rem(m,2);
        m=fix(m/2);
    end
end
if  t>=0        %处理小数点后的位数，乘2取整法 ，当乘2变为整数后结束
    for i=1:8
        b(i)=fix(t*2);
        t=2*t-fix(2*t);
    end
end
x=[a,b];
end