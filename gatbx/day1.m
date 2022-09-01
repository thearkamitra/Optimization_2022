Nind=400;           %个体数目
Maxgen=100;         %最大遗传代数
Nvar1=15;           %变量1二进制位数
Nvar2=10;           %变量2二进制位数
Nvar3=10;           %变量3二进制位数
Ggap=0.9;           %代沟
trace1=[];trace2=[];trace3=[];trace4=[];trace5=[];%性能跟踪
%建立区域描述器
Chrom=crtbp(Nind,Nvar1+Nvar2+Nvar3);       %初始种群
Field=[15,10,10;15,0.01,8;80,0.52,30;0,0,0;0,0,0;1,1,1;1,1,1];
v=bs2rv(Chrom,Field);
v1=v(:,1);
v2=v(:,2);  
v3=v(:,3);                   %初始种群十进制转换
gen=1;
while gen<Maxgen
    [Nind,N]=size(Chrom);
    M=fix(Nind/3);
    N=M;
    for i=1:N
    ObjV1t(i)=f1(v1(i),v2(i),v3(i));    
    end                                %分组后第一目标函数值
    ObjV1=ObjV1t';
    FitV1=ranking(ObjV1);              
    [m1,d] = hmdistance(FitV1);
    FitV1n=FitV1.*m1;                   %分配适应度值
    SelCh1=select('sus',Chrom(1:M,:),FitV1,Ggap);      %选择
     for i=1:N
    ObjV2t(i)=f2(v1(M+i),v2(M+i),v3(M+i));   
     end
    ObjV2=ObjV2t';                %分组后第一目标函数值
    FitV2=ranking(ObjV2);    
    [m2,d] = hmdistance(FitV2);
    FitV2n=FitV2.*m2;                   %分配适应度值
    SelCh2=select('sus',Chrom(M+1:2*M,:),FitV2,Ggap);    %选择
     for i=1:N
     ObjV3t(i)=f3(v1(2*M+i),v2(2*M+i),v3(2*M+i));  
     end                           %分组后第一目标函数值
     ObjV3=ObjV3t';   
     FitV3=ranking(ObjV3);      
     [m3,d] = hmdistance(FitV3);
    FitV3n=FitV3.*m3;                %分配适应度值
    SelCh3=select('sus',Chrom(2*M+1:3*M,:),FitV3,Ggap);   %选择
  %  for i=1:N
  %     ObjV4t(i)=f4(v1(3*M+i),v2(3*M+i)); 
  %   end                           
  %  ObjV4=ObjV4t';
  %  FitV4=ranking(ObjV4);            
  %  SelCh4=select('sus',Chrom(3*M+1:Nind,:),FitV4,Ggap); %选择
    SelCh=[SelCh1;SelCh2;SelCh3];   %合并
    Selch=recombin('xovsp',SelCh,0.7);     %重组
    Chrom=mut(SelCh);                      %变异
    v=bs2rv(Chrom,Field);        
    trace1(gen,1)=max(f1(v));
    trace1(gen,2)=sum(f1(v))/length(f1(v));
    trace2(gen,1)=max(f2(v));
    trace2(gen,2)=sum(f2(v))/length(f2(v));
    trace3(gen,1)=max(f3(v));
    trace3(gen,2)=sum(f3(v))/length(f3(v));
  %  trace4(gen,1)=min(f4(v));
  %  trace4(gen,2)=sum(f4(v))/length(f4(v));
    trace5(gen,1)=max(f1(v)+f2(v)+f3(v));
    trace5(gen,2)=sum(f1(v))/length(f1(v))+sum(f2(v))/length(f2(v))+sum(f3(v))/length(f3(v));
    gen=gen+1;
end