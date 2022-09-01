function [m,d,y] = hmdistance2(x)
[N,t]=size(x);y=zeros(N,N);d=zeros(N,1);dtemp=0;
for k=1:N   
    v(k)=x(k);
       e=tentotwo2(v(k));
for j=1:N 
    b(j)=x(j);
    c=tentotwo2(b(j));
    [~,U]=size(c);
    for i=1:U
    y(j,i)=c(1,i);
    end
        [~,W]=size(e);
        for i=1:W
            if y(j,i)~=e(1,i)
              dtemp=dtemp+1;
            end
        end
        d(j)=dtemp;
        dtemp=0;
m(k)=1/N*exp(-d(1));
if j>=2
m(k)=1/N*exp(-d(j))+m(k);
end
end
end
end

function x= tentotwo2( n )
	 m=fix(n);
	 t=n-m;  
     x=[];b=[];a=[];
	if m>=0
     for i=1:18                        % 处理整数 
	   a(19-i)=rem(m,2);
       m=fix(m/2);
     end	
    end
	if  t>0        %处理小数点后的位数，乘2取整法 ，当乘2变为整数后结束
	    for i=1:1
        b(i)=fix(t*2);
		t=2*t-fix(2*t);
        end
    end
    x=a(1:5);
end