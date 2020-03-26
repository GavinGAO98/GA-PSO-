%%%%%先以一维变量为例
%%%%%必要时可将每个环节作为函数
%%%%%可添加误差参数error，当小于误差即可跳出
%%%%%使用GA算法时，自变量区间必须离散化，本例中取步长Prcs=0.001

function [GXd,Gbest] = MyGA(fitness,Prcs,m,Pc,Pm,Cir,Xmin,Xmax)
    %%输入参数
        %fitness适应度函数
        %Prcs求解精度（步长）
        %m种群规模
        %Pc交叉概率
        %Pm变异概率
        %Cir迭代次数
        %Xmin自变量取值范围下限，是个列向量
        %Xmax自变量取值范围上限，是个列向量

    %%输出返回值
        %GXd适应度函数取最大值时最优染色体的编码
        %Gbest最大适应度值
    
    %%编码
    len = ceil(log2((Xmax-Xmin)/Prcs+1));     %计算二进制编码位数len（即码长）.根据公式：(Xmax-Xmin)/精度+1<2^len
    
    %%种群初始化
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%为使转化成二进制后数组的列数满足码长，需在前面乘以2^(n-1),且初始化一定要用ones阵
    Xd = (2^len-1)*(ones(m,1));               %Xd是m（种群规模）行十进制列向量。
    Xb = dec2bin(Xd);                       %Xb是将m个染色体的转换成二进制后的字符串矩阵，其行数是m,列数是码长，每个行向量是字符串
    XdGen = (2^len-1)*ones(m,1);              %XGen用于存储每次迭代后的染色体
    XbGen = dec2bin(XdGen);
    
    for i=1:m                               %生成m个二进制随机数储存在Xb中，使用randi()函数产生的是double型数据
        Xb(i,:) = num2str(randi([0,1],1,len),'%0d');    %将十进制转换成字符串才能赋值,格式控制符表示去掉字符串转矩阵时的空格
    end
    
    XdDcode = zeros(m,Cir+1);           %定义一个矩阵XdDcode，记录m（行）个粒子中的第i个粒子经过迭代Cir（行）次后对应的十进制自变量
    PerFit = zeros(m,Cir+1);            %定义一个矩阵PerFit，记录m（行）个粒子中的第i个粒子经过迭代Cir（行）次后的最新适应度
    
    for j=1:Cir                    %共迭代Cir次
    
        %%计算适应度值
        
        for i=1:m
            PerFit(i,j) = fitness(DecodeChrom(Xmin,Xmax,len,Xb(i,:)));
        end
    
        %%复制（选择）（目的：淘汰适应度最小的染色体）
        %%轮盘赌选择法
        pSel = zeros(m,1);                 %定义m维选择概率列向量pSel
        qAcc = zeros(m,1);                 %定义m维积累概率列向量qAcc
        FitSum = sum(PerFit(:,j));         %FitSum求种群中所有染色体适配度之和
        qTmp = 0;                          %把上一个染色体的适应度赋给暂存变量
        
        for i=1:m
            pSel(i) = PerFit(i,j)/ FitSum;
            qAcc(i) = pSel(i)+qTmp;
            qTmp = qAcc(i);                
        end

        iSel = zeros(m,1);                 %记录下轮盘赌法每次选中的（m个中）染色体编号
%         XdGen = (2^len-1)*eye(m,1);              %XdTmp用于暂存选择后的染色体
%         XbGen = dec2bin(XdGen);

        for i=1:m
            rSel = rand;                   %定义随机数
            qAddRand = [qAcc;rSel];        %将随机数添加至qAcc中组成更大的列向量qAddRand进行排序
            qAddRand = sort(qAddRand);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%向量元素升序排序后随机数所在的索引就是选择的染色体编号
            iSel(i) = find(qAddRand == rSel);
            XbGen(i,:) = Xb(iSel(i),:);    %将选择出来的适应度高的染色体赋给下一代
            
        end

        %%交叉（采用单点交叉算子）
        rPc = rand;                   %若交叉随机数rPc小于交叉概率Pc，则进行交叉
        if(rPc<Pc)
            for k=1:m/2                        %两两配对 
                rCross = randi([1,len]);       %确定交叉点位数,2个配对的染色体从第rCross位开始交叉
                   
                   %%%%%字符串截断 
                   PatLastrC = XbGen(2*k-1,rCross:len);     %提取待交叉的染色体后rCross位作为子阵      
                   MatLastrC = XbGen(2*k,rCross:len);    

                   %%%%%%%%%%%%后rCross位进行交叉，strcat(str1,st2)将str1和str2拼接上
                   XbGen(2*k-1,:) = strcat(XbGen(2*k-1,1:rCross-1),MatLastrC);
                   XbGen(2*k,:) = strcat(XbGen(2*k,1:rCross-1),PatLastrC);
            end
        end

        %%变异（采用基本位变异，对随机的一个基因或者几个基因进行变异）
        rPm = rand;                     %若交叉随机数rPc小于交叉概率Pc，则进行交叉

        if(rPm<Pm)
            for i=1:m
                rNumMut = randi([1,len]);      %从1-n中确定随机的变异基因的个数rNumMut
                rMut = randperm(len,rNumMut);  %rMut是rNumMut维行向量，记录了rNumMut个变异基因的位置
                                             %randperm(len,t)作用是生成【1，len】间t个不同的随机整数

                for t=1:rNumMut
                    %%%%%%%%%%%%%%%%%%%%%%%将随机选中的第rMut(t)位基因取反（字符型变量先用ASCII码进行整数运算后转成字符）
                    XbGen(i,rMut(t)) = num2str('1'-XbGen(i,rMut(t)));
                end
            end
        end
        
              
        %%记录每次操作后的（解码后的）自变量和适应度值，用于画图
        for i=1:m 
            Xb(i,:) = XbGen(i,:);          %将上次迭代后的染色体暂存
            XdDcode(i,j+1) = DecodeChrom(Xmin,Xmax,len,XbGen(i,:));    %第j次循环第i个染色体解码后十进制自变量XdDcode(i,j)
            PerFit(i,j+1) = fitness(XdDcode(i,j+1));       
        end
        
        %%判断是否满足精度要求，若满足可以提前跳出
%         if max(PerFit(:,j))<error
%            break; 
%         end    
    end
    
    %%循环外画图
    
    
    %%循环外解码
    Gbest = max(PerFit(:,Cir+1));     
    [GbestRow,GbestCol] = find(PerFit == Gbest,1,'last');    %找到所有适应度最大的染色体在矩阵中的位置
    GXd = XdDcode(GbestRow,GbestCol);
end

function [Xd] = DecodeChrom(Xmin,Xmax,n,Xb)
    %%函数作用：（一维变量）解码染色体，从二进制变成十进制
    %%输入参数：
        %Xmin自变量取值范围下限，是个列向量(一维时是个数)
        %Xmax自变量取值范围上限，是个列向量
        %n二进制编码位数（码长）
        %m行n列矩阵Xb,  Xb的列是一次（选择、交叉、变异）遗传操作后的染色体的二进制编码        
    %%输出参数：解码后对应于原定义域的十进制Xd
    Xd = Xmin + (Xmax-Xmin)/2^n* bin2dec(Xb);
end