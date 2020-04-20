%%先以一维为例
function [fGbest,XGPos] = MyGAPSO(fitness,Prcs,w,c1,c2,r1,r2,Pc,Pm,PopSize,D,loop,Xmin,Xmax)
    %%输入参数
        %fitness适应度函数
        %Prcs求解精度（步长）
        %w惯性权重
        %c1,c2学习因子（加速常数）
        %r1,r2（0，1）间的随机量
        %Pc交叉概率
        %Pm变异概率
        %PopSize种群规模，必须是4的倍数才可以分前后半群后再父母配对
        %D搜索空间的维数（自变量的个数） 
        %loop迭代次数
        %Xmin自变量取值范围下限，是个列向量
        %Xmax自变量取值范围上限，是个列向量
   %%输出参数
        %fGbest最后一次迭代后的全局极值
        %XGPos取得全局极值时的自变量
   
   
   %%(遗传算法初始化）编码，计算码长len
     len = ceil(log2((Xmax-Xmin)/Prcs+1));
     Xd = (2^len-1)*(zeros(PopSize/2,1));               %Xd是PopSize/2（种群规模）行十进制列向量。
     Xb = dec2bin(Xd,len);                       %Xb是将m个染色体的转换成二进制后的字符串矩阵，其行数是PopSize/2,列数是码长，每个行向量是字符串
   
   %%(用PSO的方式)初始化种群（初始化参数在目标函数m文件中完成）
   %%%%%%%%%%%%%%%%%    X是自变量 
    X = zeros(PopSize,D);       %定义一个PopSize行D列的位置矩阵，矩阵的每个行向量X(i,:)是第i个粒子的位置信息,该信息随着迭代次数而刷新
    V = zeros(PopSize,D);       %定义一个PopSize行D列的速度矩阵,
    Vmax = 0.6*Xmax;           %Vmax 通常选为 k·Xmax，其中 0.1 ≤ k ≤ 1.0，取k=0.6
    
    for i=1:PopSize
        for j=1:D
            %%% X(i,:)=Xmin+(Xmax-Xmin)*rand 对于同一个i的每个分量（变量）的rand（1）随机数是相同的，不够随机
            X(i,j) = Xmin(j)+(Xmax(j)-Xmin(j))*rand;   %%随机生成定义域内的初始位置
            V(i,j) = -Vmax(j)+2*Vmax(j)*rand;          %%随机生成定义域内的初始速度
        end
    end
    
    %%%%%%%%%%%%%%%% 定义一个数组记录第i个粒子（共PopSize个粒子）在第k次循环（共loop次循环）后的适应度值
    PerFit = zeros(PopSize,loop);
    %前一半种群PopSize/2个粒子每次迭代后个体极值（最大适应度）Indvbest，是个PopSize/2矩阵。第i个元素是第i个粒子的最大值
    Indvbest = zeros(PopSize/2,loop);      %初始化个体极值
    %PopSize/2个粒子经过loop次迭代后取得个体极值的位置IndvPos，是个PopSize/2行D列矩阵。第i行是第i个粒子的最大值的位置
    IndvPos = zeros(PopSize/2,D,loop);   %记录每次迭代后一半种群（D维变量）的位置
    
    Gbest = zeros(loop,1);      %记录每次循环后的全局极值（种群最大适应度）
    GPos = zeros(loop,D);       %记录每次循环后的取得全局极值所在的位置
    
    %%迭代loop次
    for k=1:loop
        
        %%计算适应度值
        for i=1:PopSize
            PerFit(i,k) = fitness(X(i,:));  
        end
        
        %%%%%%%%%%%满足精度要求，直接跳出循环
%         if max(PerFit(:,k)<err)
%            break;
%         end
        
        %%筛选适应度高的前一半种群
        SppFltrX = zeros(PopSize/2,D);            %记录前一半种群的位置（十进制自变量值）
        SppV = zeros(PopSize/2,D);
        [FitSort,Index] = sort(PerFit(:,k),'descend');  %将每次迭代的适应度按降序排列放在FitSort中
        
        %%更新个体极值和全局极值
        for t=1:PopSize/2
            SppFltrX(t,:) = X(Index(t),:);  %取每次适应度高的前一半种群
            SppV(t,:) = V(Index(t),:);      %取每次适应度高的前一半种群对应的速度
            
            Indvbest(t,k) = FitSort(t);     %个体极值计算
            IndvPos(t,:,k) = SppFltrX(t,:);
        
        end

        %PopSize/2个粒子的全局极值Gbest，是个常数
        %PopSize/2个粒子（经过loop次迭代后）取得全局极值的位置GPos，是个D维列向量
        Gbest(k) = FitSort(1);    
        GPos(k,:) = SppFltrX(1,:);
        
        %%用PSO公式提高前一半粒子，提高后的粒子记录在SppFltrX中
        for t=1:PopSize/2                       %%将每个粒子迭代后适应度最大的位置给
%             Vtmp = SppV(t,:);              %记录第i个粒子上一次迭代后的速度
%             Xtmp = SppFltrX(t,:);              %记录第i个粒子上一次迭代后的位置
            SppV(t,:) = w*SppV(t,:)+c1*r1*(IndvPos(t,:,k)-SppFltrX(t,:))+c2*r2*(GPos(k,:)-SppFltrX(t,:));
            SppFltrX(t,:) = SppFltrX(t,:)+SppV(t,:);
            
            %%保证每次迭代后的速度和位置在解空间内(可以有更好的算法思想来完成)
            for j=1:D
                if SppV(t,j)>Vmax(j)
                    SppV(t,j) = Vmax(j);
                elseif SppV(t,j)<(-Vmax(j))
                    SppV(t,j) = -Vmax(j);
                end
                
                if SppFltrX(t,j)>Xmax(j)
                    SppFltrX(t,j) = Xmax(j);
                elseif SppFltrX(t,j)<Xmin(j)
                    SppFltrX(t,j) = Xmin(j);
                end
            end
        
        end
        
        %%遗传操作
        
        %将前半种群编码
%         for t=1:PopSize/2                               %生成m个二进制随机数储存在Xb中，使用randi()函数产生的是double型数据
%             Xb(t,:) = CodeChrom(Prcs,Xmin,Xmax,len,SppFltrX(t,:));    %将十进制转换成字符串才能赋值,格式控制符表示去掉字符串转矩阵时的空格
%         end
        
        %随机竞争方法确定父母代
        CmptChrom = zeros(PopSize/2,D);  %随机竞争后的染色体向量
        Xtmp = SppFltrX;                    %暂存待竞争的染色体
        
        for i=1:PopSize/2-1
            rCmpt1 = unifrnd(1,PopSize/2-i+1); %四舍五入生成随机数索引
            CmptIndex1 = round(rCmpt1);
            %%%%%%%%%%%%%%%%将第一次选择出来的待比较的染色体暂存，方便最后一步添加回去
            TmpChrom1 = Xtmp(CmptIndex1,:);   
            CmptFit1 = fitness(TmpChrom1);
            Xtmp(CmptIndex1,:) = [];            %将已选择的删去避免重复竞争
            
            rCmpt2 = unifrnd(1,PopSize/2-i); %四舍五入生成随机数索引
            CmptIndex2 = round(rCmpt2);
            CmptFit2 = fitness(Xtmp(CmptIndex2,:));
            
            if(CmptFit1>CmptFit2)
               CmptChrom(i,:) = SppFltrX(CmptIndex1,:);%每次选出两个染色体比较适应度，适应度值大的作为父（母）代
            else
               CmptChrom(i,:) = SppFltrX(CmptIndex2,:);
               Xtmp(CmptIndex2,:) = [];
               Xtmp(end+1,:) = TmpChrom1;      %如果染色体2的适应度大于1，将染色体1补回
               
            end    
        end
        CmptChrom(end,:) = SppFltrX(PopSize/2,:);
        
        %将随机竞争选出的一半种群编码
        for t=1:PopSize/2                               %生成m个二进制随机数储存在Xb中，使用randi()函数产生的是double型数据
            Xb(t,:) = dec2bin((CmptChrom(t,:)-Xmin)/(Xmax-Xmin)*(2^len-1),len);    %将十进制转换成字符串才能赋值,格式控制符表示去掉字符串转矩阵时的空格
        end
        
        for n=1:PopSize/4
            PatChrom = Xb(2*n-1,:);      %适应度值大的作为父代染色体         
            MatChrom = Xb(2*n,:);
            
            %%交叉（采用单点交叉算子）
            rPc = rand;                   %若交叉随机数rPc小于交叉概率Pc，则进行交叉
            if(rPc<Pc)
                rCross = randi([1,len]);       %确定交叉点位数,2个配对的染色体从第rCross位开始交叉

                %%%%%字符串截断 
                PatLastrC = PatChrom(rCross:len);     %提取待交叉的染色体后rCross位作为子阵      
                MatLastrC = MatChrom(rCross:len);    

                %%%%%%%%%%%%后rCross位进行交叉，strcat(str1,st2)将str1和str2拼接上
                Xb(2*n-1,:) = strcat(PatChrom(1:rCross-1),MatLastrC);
                Xb(2*n,:) = strcat(MatChrom(1:rCross-1),PatLastrC);
                
            end
        end
        
        %%变异（采用基本位变异，对随机的一个基因或者几个基因进行变异）
        rPm = rand;                     %若交叉随机数rPc小于交叉概率Pc，则进行交叉

        if(rPm<Pm)
            rNumMut = randi([1,len]);      %从1-n中确定随机的变异基因的个数rNumMut
            rMut = randperm(len,rNumMut);  %rMut是rNumMut维行向量，记录了rNumMut个变异基因的位置
                                         %randperm(len,t)作用是生成【1，len】间t个不同的随机整数
            for t=1:PopSize/2
                for m=1:rNumMut
                    %%%%%%%%%%%%%%%%%%%%%%%将随机选中的第rMut(m)位基因取反（字符型变量先用ASCII码进行整数运算后转成字符）
                    Xb(t,rMut(m)) = num2str('1'-Xb(t,rMut(m)));
                end
            end
        end
        
        %%解码
%         for t=1:PopSize/2
%             Xd(t) = DecodeChrom(Xmin,Xmax,len,Xb(t,:));
%         end
        
        %%将PSO提高的前一半种群和遗传提高的后一半种群合并存放在X(PopSize,D)中
        for t=1:PopSize/2
            X(t,:) = SppFltrX(t,:);
            X(t+PopSize/2,:) = DecodeChrom(Xmin,Xmax,len,Xb(t,:));   %%后一半解码
            
        end
    end
    
    %%最后一次迭代后的适应度和位置储存在Gbest(loop+1)中
    for i=1:PopSize
            PerFit(i,loop+1) = fitness(X(i,:));  
    end
    
    Gbest(1) = [];  %将初始化时计算的最大适应度值和位置去除，最终最优解和迭代次数一一对应
    GPos(1,:) = []; 
    [Gbest(end+1),MaxIndex] = max(PerFit(:,loop+1));
    GPos(end+1,:) = X(MaxIndex,:);
    
    fGbest = Gbest(loop);
    XGPos = GPos(loop,:);
    
    %%作图描述
    %作出最佳适应度值和迭代次数的曲线
   
        plot(Gbest,'-*');
        xlabel('Generations');
        ylabel('Best fitness');
        title('最佳适应度值关于迭代次数的变化图');
        
    %以迭代次数为参数（即图上有loop个点）作出最优解位置变化（以各个分量为坐标系）曲线
%         Xaxis = linspace(1,loop,loop);
%         plot(Xaxis,GPos);
%         
    %作出任务分配各点的路径
    
end

function [Xd] = DecodeChrom(Xmin,Xmax,len,Xb)
    %%函数作用：（一维变量）解码染色体，从二进制变成十进制
    %%输入参数：
        %Xmin自变量取值范围下限，是个列向量(一维时是个数)
        %Xmax自变量取值范围上限，是个列向量
        %len二进制编码位数（码长）
        %PopSize行len列矩阵Xb,  Xb的列是一次（选择、交叉、变异）遗传操作后的染色体的二进制编码        
    %%输出参数：Xd解码后对应于原定义域的十进制
    Xd = Xmin + (Xmax-Xmin)/(2^len-1)* bin2dec(Xb);
end
