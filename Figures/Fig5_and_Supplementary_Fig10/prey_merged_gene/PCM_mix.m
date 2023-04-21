clear all;

%% load data

%type = 'full';
%type = 'half';
type = 'quarter';

filename = ['data_merge_',type];
load(filename)

lag = 25;
maxL= length(y(:,1));
E = 3; 
tau = lag;
N = 4; 
G = y;

LL=zeros(E,maxL,N);

for i=1:N
    for p=1+(E-1)*tau:maxL
        for q=1:E
            LL(q,p,i)=G(p-(q-1)*tau,i);
        end
    end
end

MM=zeros(maxL,maxL,N);
for i=1:N
    for r=1+(E-1)*tau:maxL
        for s=r:maxL
            MM(r,s,i)=norm(LL(:,r,i)-LL(:,s,i));
        end
    end
    MM(:,:,i)=MM(:,:,i)+MM(:,:,i)';
end
GG=G(1+(E-1)*tau:maxL,:);
MM=MM(1+(E-1)*tau:maxL,1+(E-1)*tau:maxL,:);

%被预测序列的长度
%PL=5;
PL=maxL - (E-1)*tau - lag*2-1;
%搜寻时滞的大小
maxfdelay= lag;

fdelay=0:1:maxfdelay;
%被预测序列的指标
cy=maxfdelay+1:maxfdelay+PL;
%Library的大小
LB=PL;
%忽略时间近邻的个数
NB=0;

C=zeros(length(fdelay),N,N);
CX=zeros(PL,length(fdelay),N,N);
CY=zeros(PL,N,N);
PC=zeros(length(fdelay),N,N,N);
PCX=zeros(PL,length(fdelay),N,N,N);
PCY=zeros(PL,N,N,N);
PPC=zeros(N,N);
PV=zeros(N,N);


%%%%% Gi预测Gj %%%%%
for i=1:N
    for j=1:N
        if j==i
            continue;
        end
        %被预测序列
        CY(:,i,j)=GG(cy,j);
        sh=floor(rand(1,1)*(maxL-(E-1)*tau-maxfdelay-(maxfdelay+LB)));
        for w=1:length(fdelay)
            %Library的指标
            mx=maxfdelay+1+sh:maxfdelay+LB+sh;
            for r=1:PL
                %用来预测的点在原序列中的指标
                px=maxfdelay+r+fdelay(w);
                %Library到此点的距离序列
                MX=MM(mx,px,i);
                %找寻最近且非时间邻近的Library点，其在原序列中的指标存入IX中
                IX=zeros(E+1,1);
                [~,ind]=sort(MX,'ascend');
                for t=1:LB
                    if MX(ind(t))~=0
                        break;
                    end
                end
                tt=0;
                fd=0;
                while tt<E+1
                    if abs(ind(t+fd)+maxfdelay+sh-px)<=NB
                        fd=fd+1;
                    else
                        IX(tt+1)=ind(t+fd)+maxfdelay+sh;
                        tt=tt+1;
                        fd=fd+1;
                    end
                end
                %计算权重
                U=zeros(E+1,1);
                sumU=0;
                for b=1:E+1
                    U(b)=exp(-MX(IX(b)-maxfdelay-sh)/MX(IX(1)-maxfdelay-sh));
                    sumU=sumU+U(b);
                end
                %计算预测点值
                sumY=0;
                for b=1:E+1
                    sumY=sumY+(U(b)/sumU)*GG(IX(b)-fdelay(w),j);
                end
                CX(r,w,i,j)=sumY;
            end
            CR=corrcoef(CX(:,w,i,j),CY(:,i,j));
            C(w,i,j)=CR(1,2);
        end
    end
end

%%%%% PGj_Gi预测G_k %%%%%
for i=1:N
    for j=1:N
        if i==j
            continue;
        end
        [mv,mid]=max(C(:,i,j));
        P=CX(:,mid,i,j);
        PP=[zeros((E-1)*tau+maxfdelay,1);P;zeros(length(GG(:,1))-maxfdelay-PL,1)];
        LLP=zeros(E,maxL);
        for p=1+(E-1)*tau:maxL
            for q=1:E
                LLP(q,p)=PP(p-(q-1)*tau);
            end
        end
        MMP=zeros(maxL,maxL);
        for r=1+(E-1)*tau:maxL
            for s=r:maxL
                MMP(r,s)=norm(LLP(:,r)-LLP(:,s));
            end
        end
        MMP=MMP+MMP';
        PP=PP(1+(E-1)*tau:maxL);
        MMP=MMP(1+(E-1)*tau:maxL,1+(E-1)*tau:maxL);
        for k=1:N
            if k==i || k==j
                continue;
            end
            %被预测序列
            PCY(:,k,i,j)=GG(cy,k);
            sh=floor(rand(1,1)*(maxL-(E-1)*tau-maxfdelay-(maxfdelay+LB)));
            for w=1:length(fdelay)
                %Library的指标
                mx=maxfdelay+1+sh:maxfdelay+LB+sh;
                for r=1:PL
                    %用来预测的点在原序列中的指标
                    px=maxfdelay+r+fdelay(w);
                    %Library到此点的距离序列
                    MX=MMP(mx,px);
                    %找寻最近且非时间邻近的Library点，其在原序列中的指标存入IX中
                    IX=zeros(E+1,1);
                    [~,ind]=sort(MX,'ascend');
                    for t=1:LB
                        if MX(ind(t))~=0
                            break;
                        end
                    end
                    tt=0;
                    fd=0;
                    while tt<E+1
                        if abs(ind(t+fd)+maxfdelay+sh-px)<=NB
                            fd=fd+1;
                        else
                            IX(tt+1)=ind(t+fd)+maxfdelay+sh;
                            tt=tt+1;
                            fd=fd+1;
                        end
                    end
                    %计算权重
                    U=zeros(E+1,1);
                    sumU=0;
                    for b=1:E+1
                        U(b)=exp(-MX(IX(b)-maxfdelay-sh)/MX(IX(1)-maxfdelay-sh));
                        sumU=sumU+U(b);
                    end
                    %计算预测点值
                    sumY=0;
                    for b=1:E+1
                        sumY=sumY+(U(b)/sumU)*GG(IX(b)-fdelay(w),k);
                    end
                    PCX(r,w,k,i,j)=sumY;
                end
                CR=corrcoef(PCX(:,w,k,i,j),PCY(:,k,i,j));
                PC(w,k,i,j)=CR(1,2);
            end
        end
    end
end

for i=1:N
    for j=1:N
        if j==i
            continue;
        end
        PA=[];
        for k=1:N
            if k==i || k==j
                continue;
            end
            [mv,mid]=max(PC(:,j,i,k));
            PPA=PCX(:,mid,j,i,k);
            PA=[PA PPA];
        end
        [mv,mid]=max(C(:,i,j));
        [PPC(i,j),PV(i,j)]=partialcorr(CX(:,mid,i,j),CY(:,i,j),PA);
    end
end

filename = ['PCM_result_',type];
save(filename, 'PV','PPC')