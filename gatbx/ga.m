clc
clear
Num=400;            %个体数目
Maxgen=50;         %最大遗传代数
Nvar1=15;           %变量1二进制位数
Nvar2=15;           %变量2二进制位数
Nvar3=10;           %变量3二进制位数
Nvar4=10;           %变量4二进制位数
Nvar5=10;           %变量5二进制位数
Nvar6=10;           %变量6二进制位数
Ggap=0.9;          %代沟
% trace1=[];trace2=[];trace3=[];trace4=[];trace5=[]; 性能跟踪
%建立区域描述器
Chrom=crtbp(Num,Nvar1+Nvar2+Nvar3+Nvar4+Nvar5+Nvar6);       %初始种群
Field=[15,15,10,10,10,10;15,15,0.01,0.01,0.01,8;50,50,0.12,0.12,0.12,30;0,0,0,0,0,0;0,0,0,0,0,0;1,1,1,1,1,1;1,1,1,1,1,1];%Em1,Em2,L1,L2,L3,Ek0
v=bs2rv(Chrom,Field);                %初始种群十进制转换
gen=1;TheBest1=[];TheBest2=[];TheBest3=[];    %最优值初始化
while gen<Maxgen
    [Nind,~]=size(Chrom);
    M=fix(Nind/3);
    if M==0
        break;
    end
    N=M;               %根据目标函数个数均分初始种群为子种群
    [Z, Xcluster1, A, cluster]= wfIsodata_ND(v(1:N,:),5,2,40,10,5,0.4,0,100); %子种群聚类
    [~,wow1]=size(Xcluster1); ObjV14n=[];ObjV11n=[];ObjV12n=[];ObjV13n=[];
    ObjV11=[];ObjV12=[];ObjV13=[];ObjV14=[];m11=[];m11t=[];m12=[];m12t=[];m13=[];m13t=[];
    if gen>1
        Xcluster1{1,1}(1,:)=TheBest1;
    end
    %%%%%%%%%%计算目标函数，每个类分配适应度%%%%%%%%%%
    if wow1>=1
        [ooh1,~]=size(Xcluster1{1,1});
        for i=1:ooh1
            ObjV1t1(i)=f0(Xcluster1{1,1}(i,1),Xcluster1{1,1}(i,2),Xcluster1{1,1}(i,3),Xcluster1{1,1}(i,4),Xcluster1{1,1}(i,5),Xcluster1{1,1}(i,6));
        end
        ObjV11=ObjV1t1';
        [m11] = hmdistance1(ObjV11); m11t=m11';
        ObjV11n=ObjV11.*m11t;
    end
    if wow1>=2
        [pph1,~]=size(Xcluster1{1,2});
        for i=1:pph1
            ObjV1t2(i)=f0(Xcluster1{1,2}(i,1),Xcluster1{1,2}(i,2),Xcluster1{1,2}(i,3),Xcluster1{1,2}(i,4),Xcluster1{1,2}(i,5),Xcluster1{1,2}(i,6));
        end
        ObjV12=ObjV1t2';
        [m12] = hmdistance1(ObjV12); m12t=m12';
        ObjV12n=ObjV12.*m12t;
    end
    if wow1>=3
        [qqh1,~]=size(Xcluster1{1,3});
        for i=1:qqh1
            ObjV1t3(i)=f0(Xcluster1{1,3}(i,1),Xcluster1{1,3}(i,2),Xcluster1{1,3}(i,3),Xcluster1{1,3}(i,4),Xcluster1{1,3}(i,5),Xcluster1{1,3}(i,6));
        end
        ObjV13=ObjV1t3';
        [m13] = hmdistance1(ObjV13); m13t=m13';
        ObjV13n=ObjV13.*m13t;
    end
    if wow1>=4
        [rrh1,~]=size(Xcluster1{1,4});
        for i=1:rrh1
            ObjV1t4(i)=f0(Xcluster1{1,4}(i,1),Xcluster1{1,4}(i,2),Xcluster1{1,4}(i,3),Xcluster1{1,4}(i,4),Xcluster1{1,4}(i,5),Xcluster1{1,4}(i,6));
        end
        ObjV14=ObjV1t4';
        [m14] = hmdistance1(ObjV14); m14t=m14';
        ObjV14n=ObjV14.*m14t;
    end
    ObjV1=[ObjV11n;ObjV12n;ObjV13n;ObjV14n];
    Chrom1=tr(Xcluster1);                             %将十进制转换为二进制
    %%%%%%保留最优值%%%%%%%
    [~,i1] = max(ObjV1);
    if   i1<=ooh1
        TheBest1=Xcluster1{1,1}(i1,:);
    end
    if   i1>ooh1&&i1<=(ooh1+pph1)
        b=i1-ooh1;
        TheBest1=Xcluster1{1,2}(b,:);
    end
    if   i1>(ooh1+pph1)&&i1<=(qqh1+ooh1+pph1)
        b=i1-ooh1-pph1;
        TheBest1=Xcluster1{1,3}(b,:);
    end
    SelCh1=select('sus',Chrom1,ObjV1,Ggap);           %选择
    [Z, Xcluster2, A, cluster]= wfIsodata_ND(v(N+1:2*N,:),5,2,40,10,5,0.4,0,100);   %子种群聚类
    [~,wow2]=size(Xcluster2); ObjV21n=[];ObjV23n=[];ObjV22n=[];ObjV24n=[];
    ObjV21=[];ObjV22=[];ObjV23=[];ObjV24=[];m21=[];m22=[];m23=[];m21t=[];m22t=[];m23t=[];
    if gen>1
        Xcluster2{1,1}(1,:)=TheBest2;
    end
    %%%%%%%%%%计算目标函数，每个类分配适应度%%%%%%%%%%
    if wow2>=1
        [ooh2,~]=size(Xcluster2{1,1});
        for i=1:ooh2
            ObjV2t1(i)=f0(Xcluster2{1,1}(i,1),Xcluster2{1,1}(i,2),Xcluster2{1,1}(i,3),Xcluster2{1,1}(i,4),Xcluster2{1,1}(i,5),Xcluster2{1,1}(i,6));
        end
        ObjV21=ObjV2t1';
        [m21] = hmdistance2(ObjV21); m21t=m21';
        ObjV21n=ObjV21.*m21t;
    end
    if wow2>=2
        [pph2,~]=size(Xcluster2{1,2});
        for i=1:pph2
            ObjV2t2(i)=f0(Xcluster2{1,2}(i,1),Xcluster2{1,2}(i,2),Xcluster2{1,2}(i,3),Xcluster2{1,2}(i,4),Xcluster2{1,2}(i,5),Xcluster2{1,2}(i,6));
        end
        ObjV22=ObjV2t2';
        [m22] = hmdistance2(ObjV22); m22t=m22';
        ObjV22n=ObjV22.*m22t;
    end
    if wow2>=3
        [qqh2,~]=size(Xcluster2{1,3});
        for i=1:qqh2
            ObjV2t3(i)=f0(Xcluster2{1,3}(i,1),Xcluster2{1,3}(i,2),Xcluster2{1,3}(i,3),Xcluster2{1,3}(i,4),Xcluster2{1,3}(i,5),Xcluster2{1,3}(i,6));
        end
        ObjV23=ObjV2t3';
        [m23] = hmdistance2(ObjV23); m23t=m23';
        ObjV23n=ObjV23.*m23t;
    end
    if wow2>=4
        [rrh2,~]=size(Xcluster2{1,4});
        for i=1:rrh2
            ObjV2t4(i)=f0(Xcluster2{1,3}(i,1),Xcluster2{1,3}(i,2),Xcluster2{1,3}(i,3),Xcluster2{1,3}(i,4),Xcluster2{1,3}(i,5),Xcluster2{1,3}(i,6));
        end
        ObjV24=ObjV2t4';
        [m24] = hmdistance2(ObjV24); m24t=m24';
        ObjV24n=ObjV24.*m24t;
    end
    ObjV2d=[ObjV21n;ObjV22n;ObjV23n;ObjV24n];
    ObjV2x=ds(ObjV2d)*0.0000001;         %转化为最小值优化问题
    %%%%%%保留最优值%%%%%%%
    [~,i2] = max(ObjV2x);
    if   i2<=ooh2
        TheBest2=Xcluster2{1,1}(i2,:);
    end
    if   i2>ooh2&&i2<=(pph2+ooh2)
        b=i2-ooh2;
        TheBest2=Xcluster2{1,2}(b,:);
    end
    if   i2>(pph2+ooh2)&&i2<=(qqh2+ooh2+pph2)
        b=i2-pph2-ooh2;
        TheBest2=Xcluster2{1,3}(b,:);
    end
    Chrom2=tr(Xcluster2);                            %将十进制转换为二进制
    SelCh2=select('sus',Chrom2,ObjV2x,Ggap);    %选择
    [Z, Xcluster3, A, cluster]= wfIsodata_ND(v(N+1:2*N,:),5,2,40,10,5,0.4,0,100); %子种群聚类
    [~,wow3]=size(Xcluster3); ObjV31n=[];ObjV32n=[];ObjV33n=[];ObjV34n=[];
    ObjV31=[];ObjV32=[];ObjV33=[];ObjV34=[];m31=[];m32=[];m33=[];m31t=[];m32t=[];m33t=[];
    if gen>1
        Xcluster3{1,1}(1,:)=TheBest3;
    end
    %%%%%%%%%%计算目标函数，每个类分配适应度%%%%%%%%%%
    if wow3>=1
        [ooh3,~]=size(Xcluster3{1,1});
        for i=1:ooh3
            ObjV3t1(i)=f0(Xcluster3{1,1}(i,1),Xcluster3{1,1}(i,2),Xcluster3{1,1}(i,3),Xcluster3{1,1}(i,4),Xcluster3{1,1}(i,5),Xcluster3{1,1}(i,6));
        end
        ObjV31=ObjV3t1';
        [m31] = hmdistance3(ObjV31); m31t=m31';
        ObjV31n=ObjV31.*m31t;
    end
    if wow3>=2
        [pph3,~]=size(Xcluster3{1,2});
        for i=1:pph3
            ObjV3t2(i)=f0(Xcluster3{1,2}(i,1),Xcluster3{1,2}(i,2),Xcluster3{1,2}(i,3),Xcluster3{1,2}(i,4),Xcluster3{1,2}(i,5),Xcluster3{1,2}(i,6));
        end
        ObjV32=ObjV3t2';
        [m32] = hmdistance3(ObjV32); m32t=m32';
        ObjV32n=ObjV32.*m32t;
    end
    if wow3>=3
        [qqh3,~]=size(Xcluster2{1,3});
        for i=1:qqh3
            ObjV3t3(i)=f0(Xcluster3{1,3}(i,1),Xcluster3{1,3}(i,2),Xcluster3{1,3}(i,3),Xcluster3{1,3}(i,4),Xcluster3{1,3}(i,5),Xcluster3{1,3}(i,6));
        end
        ObjV33=ObjV3t3';
        [m33] = hmdistance3(ObjV33); m33t=m33';
        ObjV33n=ObjV33.*m33t;
    end
    if wow3>=4
        [rrh3,~]=size(Xcluster3{1,4});
        for i=1:rrh3
            ObjV3t4(i)=f0(Xcluster3{1,4}(i,1),Xcluster3{1,4}(i,2),Xcluster3{1,4}(i,4),Xcluster3{1,4}(i,4),Xcluster3{1,4}(i,5),Xcluster3{1,4}(i,6));
        end
        ObjV34=ObjV3t4';
        [m34] = hmdistance3(ObjV34); m34t=m34';
        ObjV34n=ObjV34.*m34t;
    end
    ObjV3=[ObjV31n;ObjV32n;ObjV33n;ObjV34n];
    Chrom3=tr(Xcluster3);                     %将十进制转换为二进制
    ObjV3x=ds(ObjV3);          %转化为最小值优化问题
    %%%%%%保留最优值%%%%%%%
    [~,i3] = max(ObjV3x);
    if   i3<=ooh3
        TheBest3=Xcluster3{1,1}(i3,:);
    end
    if   i3>ooh3&&i3<=(ooh3+pph3)
        b=i3-ooh3;
        TheBest3=Xcluster3{1,2}(b,:);
    end
    if   i3>(pph3+ooh3)&&i3<=(ooh3+pph3+qqh3)
        b=i3-pph3-ooh3;
        TheBest3=Xcluster3{1,3}(b,:);
    end
    SelCh3=select('sus',Chrom3, ObjV3x,Ggap);   %选择
    SelCh=[SelCh1;SelCh2;SelCh3];          %合并
    SelCh=recombin('xovsp',SelCh,0.7);     %重组
    Chrom=mut(SelCh);                      %变异
    if gen<(Maxgen-1)
        v=[];ObjV1t1=[];ObjV1t2=[];ObjV1t3=[];ObjV1t4=[];
        ObjV2t1=[];ObjV2t2=[];ObjV2t3=[];ObjV2t4=[];
        ObjV3t1=[];ObjV3t2=[];ObjV3t3=[];ObjV3t4=[];
    end
    v=bs2rv(Chrom,Field);                  %原始数组更新
    %     trace1(gen,1)=max(f1(v1(i),v2(i),v3(i)));
    %     trace1(gen,2)=sum(f1(v1(i),v2(i),v3(i)))/length(f1(v1(i),v2(i),v3(i)));
    %     trace2(gen,1)=max(f2(v1(M+i),v2(M+i),v3(M+i)));
    %     trace2(gen,2)=sum(f2(v1(M+i),v2(M+i),v3(M+i)))/length(f2(v1(M+i),v2(M+i),v3(M+i)));
    %     trace3(gen,1)=max(f3(v1(2*M+i),v2(2*M+i),v3(2*M+i)));
    %     trace3(gen,2)=sum(f3(v1(2*M+i),v2(2*M+i),v3(2*M+i)))/length(f3(v1(2*M+i),v2(2*M+i),v3(2*M+i)));
    %     trace4(gen,1)=min(f4(v));
    %     trace4(gen,2)=sum(f4(v))/length(f4(v));
    %     trace5(gen,1)=max(f1(v1(i),v2(i),v3(i))+f2(v1(M+i),v2(M+i),v3(M+i)))+f3(v1(2*M+i),v2(2*M+i),v3(2*M+i));
    %     trace5(gen,2)=sum(f1(v1(i),v2(i),v3(i)))/length(f1(v1(i),v2(i),v3(i)))+sum(f2(v1(M+i),v2(M+i),v3(M+i)))/length(f2(v1(M+i),v2(M+i),v3(M+i)))+sum(f3(v1(2*M+i),v2(2*M+i),v3(2*M+i)))/length(f3(v1(2*M+i),v2(2*M+i),v3(2*M+i)));
    gen=gen+1;
end