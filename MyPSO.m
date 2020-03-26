%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%调用前必读
%切记修改搜索空间维数D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%待优化问题
%（1）求点作图、求最小值函数、部分算法改进
% (2)可以把每个环节作为函数单独提炼

%定义PSO函数(function必须是第一行)        
function [GPos,Gbest] = MyPSO(fitness,w,c1,c2,r1,r2,m,D,Cir,Xmin,Xmax)
    %%输入参数
        %适应度函数fitness
        %惯性权重w
        %学习因子(加速常数)c1,c2
        %（0，1）间的随机量r1,r2
        %粒子群规模m
        %搜索空间的维数D（即变量的个数）
        %迭代次数Cir
        %Xmin和Xmax是目标函数定义域的下（上）限D维行（或列）向量
    %%返回值
        %最终最优解的位置Xbest（是个D维的列向量）
        %最终最优解的适应度Fbest
    
    %%初始化种群（初始化参数在目标函数m文件中完成）
    X = zeros(m,D);       %定义一个m行D列的位置矩阵，矩阵的每个行向量是第i个粒子的位置信息,该信息随着迭代次数而刷新
    V = zeros(m,D);       %定义一个m行D列的速度矩阵,
    Vmax = 0.6*Xmax;      %Vmax 通常选为 k·Xmax，其中 0.1 ≤ k ≤ 1.0，取k=0.6
    
    for i=1:m
        for j=1:D
            %%% X(i,:)=Xmin+(Xmax-Xmin)*rand 对于同一个i的每个分量（变量）的rand（1）随机数是相同的，不够随机
            X(i,j) = Xmin(j)+(Xmax(j)-Xmin(j))*rand;   %%随机生成定义域内的初始位置
            V(i,j) = -Vmax(j)+2*Vmax(j)*rand;          %%随机生成定义域内的初始速度
        end
    end
    
    %%计算初始适应度
    PerFit = zeros(m,1);    %定义一个列向量PerFit，记录m（行）个粒子中的第i个粒子经过迭代Cir（行）次后的最新适应度
    
    for i=1:m
        PerFit(i) = fitness(X(i,:));   %计算每个粒子的适应度,  初始计算算第1次迭代
    end
    
    %%计算初始极值和位置
    
    %m个粒子经过Cir次迭代后个体极值（最大适应度）Indvbest，是个m维列向量。第i个元素是第i个粒子的最大值
    %m个粒子经过Cir次迭代后取得个体极值的位置IndvPos，是个m行D列矩阵。第i行是第i个粒子的最大值的位置
    %m个粒子经过Cir次迭代后全局极值Gbest，是个常数
    %m个粒子经过Cir次迭代后取得全局极值的位置GPos，是个D维列向量
    Indvbest = PerFit;      %初始化个体极值
    IndvPos = X;
    Gbest = max(Indvbest);    
    GPos = IndvPos(find(Indvbest == Gbest),:);   %find(数组==max（数组）)返回数组中最大值索引index
    
    %%第i个粒子迭代Cir次,更新速度和位置
    for k=2:Cir+1       
        for i=1:m                       %%将每个粒子迭代后适应度最大的位置给
            Vtmp = V(i,:);              %记录第i个粒子上一次迭代后的速度
            Xtmp = X(i,:);              %记录第i个粒子上一次迭代后的位置
            V(i,:) = w*V(i,:)+c1*r1*(IndvPos(i,:)-X(i,:))+c2*r2*(GPos-X(i,:));
            X(i,:) = X(i,:)+V(i,:);
            
            %%保证每次迭代后的速度和位置在解空间内(可以有更好的算法思想来完成)
            for j=1:D
                if V(i,j)>Vmax(j)
                    V(i,j) = Vmax(j);
                elseif V(i,j)<(-Vmax(j))
                    V(i,j) = -Vmax(j);
                end
                
                if X(i,j)>Xmax(j)
                    X(i,j) = Xmax(j);
                elseif X(i,j)<Xmin(j)
                    X(i,j) = Xmin(j);
                end
            end
            
            PerFit(i) = fitness(X(i,:));
            
            if PerFit(i)>Indvbest(i)    %%如果第i个粒子在第k次迭代的个体适应度大于其个体极值，则位置和适应度更新
                Indvbest(i) = PerFit(i);
                IndvPos(i,:) = X(i,:);
            end
            
            if PerFit(i)>Gbest         %%如果第i个粒子在第k次迭代的个体适应度大于全局极值，则位置和适应度更新
                Gbest = PerFit(i);
                GPos = X(i,:);
            end
           %%可以添加误差判断跳出 
        end
    end
    
end