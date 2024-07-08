
%---------------------------------------------------------------------------------------------------------------------------

 % Hunger Games Search (HGS)
% Visions, Conception, Implementation, Deep Analysis, Perspectives, and Towards Performance Shifts
% Website of HGS: http://www.aliasgharheidari.com/HGS.html

%  Yutao Yang and Professor Huiling Chen (citations above 7000) and Ali Asghar Heidari (citations above 4000)
%  College of Computer Science and Artificial Intelligence, Wenzhou University, Wenzhou, Zhejiang 325035, China
%  Exceptionally Talented Researcher, Department of Computer Science, School of Computing, National University of Singapore, Singapore
%  National Elite, Exceptionally Talented Researcher, School of Surveying and Geospatial Engineering, College of Engineering, University of Tehran, Tehran 1439957131, Iran

%  Thanks to Professor Amir H Gandomi (citations above 19000)
%  Faculty of Engineering & Information Technology, University of Technology Sydney, NSW 2007, Australia

%  Last update: 02-17-2021

%  e-Mail: yutaoyang.style@foxmail.com,
%  e-Mail: as_heidari@ut.ac.ir, aliasghar68@gmail.com,
%  e-Mail (Singapore): aliasgha@comp.nus.edu.sg, t0917038@u.nus.edu
%---------------------------------------------------------------------------------------------------------------------------
%  Co-author: Yutao Yang(yutaoyang.style@foxmail.com), Huiling Chen(chenhuiling.jlu@gmail.com), Ali Asghar Heidari(as_heidari@ut.ac.ir), Amir H Gandomi(a.h.gandomi@gmail.com)
%---------------------------------------------------------------------------------------------------------------------------

% After use, please refer to the main paper:
% Yutao Yang, Huiling Chen, Ali Asghar Heidari, Amir H Gandomi, 
% Hunger Games Search: Visions, Conception, Implementation, Deep Analysis, Perspectives, and Towards Performance Shifts
% Expert Systems With Applications, https://doi.org/10.1016/j.eswa.2021.114864 (Q1, 5-Year Impact Factor: 5.448, H-INDEX: 184)
%---------------------------------------------------------------------------------------------------------------------------

%  Researchgate: https://www.researchgate.net/profile/Ali_Asghar_Heidari
%  Website of HGS: http://www.aliasgharheidari.com/HGS.html

% You can also use and compare with our other new optimization methods: (HHO)-2019- http://www.aliasgharheidari.com/HHO.html
%                                                                       (SMA)-2020- http://www.aliasgharheidari.com/SMA.html
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Destination_fitness,bestPositions,Convergence_curve,position_history]=RMHGS(fhd,dim,N,Max_iter,lb,ub,varargin)
disp('HGS is now tackling your problem')

tic
% initialize position
bestPositions=zeros(1,dim);
tempPosition=zeros(N,dim);

Destination_fitness=inf;%change this to -inf for maximization problems
Worstest_fitness=-inf;
AllFitness = inf*ones(N,1);%record the fitness of all positions
VC1 = ones(N,1);%record the variation control of all positions

lb=lb.*ones(1,dim);
ub=ub.*ones(1,dim);

weight3 = ones(N,dim);%hungry weight of each position
weight4 = ones(N,dim);%hungry weight of each position

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);
Convergence_curve=zeros(1,Max_iter);
it=1; %Number of iterations

hungry = zeros(1,size(X,1));%record the hungry of all positions
count=0;

BestX=zeros(N,dim);
Bestfit=inf*ones(N,1);
position_history=zeros(N,Max_iter,dim);
% Main loop
while  it <= Max_iter
    VC2 = 0.03; %The variable of variation control 

    sumHungry = 0;%record the sum of each hungry 
    
    M1=mode(BestX);
    MM1=median(BestX);
    %sort the fitness
    for i=1:size(X,1)
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        mean1=(BestX(i,:)+M1+MM1)/3;
        mean1=(bestPositions+mean1)/2;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+mean1.*Flag4ub+mean1.*Flag4lb;
        position_history(i,it,:)=X(i,:);
        AllFitness(i) = feval(fhd,X(i,:)',varargin{:});
        if AllFitness(i) < Bestfit(i)
            BestX(i,:)=X(i,:);
            Bestfit(i)=AllFitness(i);
        end
    end
    [AllFitnessSorted,IndexSorted] = sort(AllFitness);
    bestFitness = AllFitnessSorted(1);
    worstFitness = AllFitnessSorted(size(X,1));
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(IndexSorted(1),:);
        Destination_fitness = bestFitness;
        count=0;
    end
    
    if worstFitness > Worstest_fitness
        Worstest_fitness = worstFitness;
    end
    
    X_DO=bestPositions;
    r1=rand();
    r2=rand();
    M=mode(BestX);
    MM=median(BestX);
    MMM=mean(BestX);
    mean2=(M+MM+MMM)/3;
    mean2=(bestPositions+mean2)/2;
    for j=1:dim
        X_DO(j)=bestPositions(j)+r1*(r2*(lb(j)+ub(j)-bestPositions(j))-bestPositions(j));
        if X_DO(j)<lb(j) || X_DO(j)>ub(j) 
            X_DO(j)=mean2(j);
        end
    end
    newfit=feval(fhd,X_DO',varargin{:});
    if newfit<Destination_fitness
        bestPositions=X_DO;
        Destination_fitness = newfit;
        count=1;
        tempPosition(count,:)=X_DO;
    end
    
    MM=ceil(N*(0.5));
    aaa=1-sin((pi/2)*(it/Max_iter));
    for i=1:MM
        A= randi([1,ceil(N*(0.5))]);
        B= randi([1,ceil(N*(0.5))]);
        if rand<aaa
            [X1,X2]=UniformCrossover(X(IndexSorted(A),:),X(IndexSorted(B),:));
            if rand>0.5
                X(IndexSorted(N+1-i),:)=X1;
            else
                X(IndexSorted(N+1-i),:)=X2;
            end
            pp=0.5*cos((pi/2)*(it/Max_iter));
            if rand<pp
                a1=rand();
                a2=rand();
                for j=1:dim
                    X(IndexSorted(N+1-i),j)=X(IndexSorted(N+1-i),j)+a1*(a2*(lb(j)+(rand*(ub(j)-lb(j))))-X(IndexSorted(N+1-i),j));
                    if X(IndexSorted(N+1-i),j)<lb(j) || X(IndexSorted(N+1-i),j)>ub(j) 
                        X(IndexSorted(N+1-i),j)=lb(j)+rand*(ub(j)-lb(j));
                    end
                end
            end
        else
            [X1,X2]=UniformCrossover(bestPositions,X(IndexSorted(A),:));
            if rand>0.5
                X(IndexSorted(N+1-i),:)=X1;
            else
                X(IndexSorted(N+1-i),:)=X2;
            end
            pp=0.5*cos((pi/2)*(it/Max_iter));
            if rand<pp
                a1=rand();
                a2=rand();
                for j=1:dim
                    X(IndexSorted(N+1-i),j)=X(IndexSorted(N+1-i),j)+a1*(a2*(lb(j)+(rand*(ub(j)-lb(j))))-X(IndexSorted(N+1-i),j));
                    if X(IndexSorted(N+1-i),j)<lb(j) || X(IndexSorted(N+1-i),j)>ub(j) 
                        X(IndexSorted(N+1-i),j)=lb(j)+rand*(ub(j)-lb(j));
                    end
                end
            end
        end  
    end
    
    for i = 1:size(X,1)
         %calculate the variation control of all positions
         VC1(i) = sech(abs(AllFitness(i)-Destination_fitness));    
         %calculate the hungry of each position
        if Destination_fitness == AllFitness(i)
            hungry(1,i) = 0;
            count = count+1;
            tempPosition(count,:)=X(i,:);
        else
            temprand = rand();
            c = (AllFitness(i)-Destination_fitness)/(Worstest_fitness-Destination_fitness)*temprand*2*(ub-lb);
            if c<100
                b=100*(1+temprand);
            else
                b=c;
            end   
            hungry(1,i) = hungry(1,i)+ max(b); 
            sumHungry = sumHungry + hungry(1,i);
        end
    end 
    
    %calculate the hungry weight of each position
    for i=1:size(X,1)
        for j=2:size(X,2)
                weight3(i,j) = (1-exp(-abs(hungry(1,i)-sumHungry)))*rand()*2;
                if rand()<VC2
                    weight4(i,j) = hungry(1,i)*size(X,1)/sumHungry*rand();
                else
                    weight4(i,j) = 1;
                end
        end
        
    end
    
    
    % Update the Position of search agents
    if it<(1/3)*Max_iter
        shrink=2*(1-it/Max_iter); % a decreases linearly fron 2 to 0
        for i=1:size(X,1)
            if  rand<VC2
                X(i,:) = X(i,j)*(1+randn(1));
            else
                A = randi([1,count]);
                B = randi([1,N]);
                C = randi([1,N]);
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    X(i,:) = weight4(i,:).*BestX(i,:)+vb.*weight3(i,:).*abs(2*tempPosition(A,:)-(X(B,:)+X(C,:)));
                else
                    X(i,:) = weight4(i,:).*BestX(i,:)-vb.*weight3(i,:).*abs(2*tempPosition(A,:)-(X(B,:)+X(C,:)));
                end
            end
        end
    elseif  it>=(1/3)*Max_iter && it<=(2/3)*Max_iter 
        shrink=2*(1-it/Max_iter); % a decreases linearly fron 2 to 0
        for i=((N/2)+1):N
            if  rand<VC2
                X(IndexSorted(i),:) = X(IndexSorted(i),j)*(1+randn(1));
            else
                A = randi([1,count]);
                B = randi([1,N*0.5]);
                C = randi([1,N*0.5]);
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    X(IndexSorted(i),:) = weight4(IndexSorted(i),:).*BestX(IndexSorted(i),:)+vb.*weight3(IndexSorted(i),:).*abs(tempPosition(A,:)-((X(IndexSorted(B),:)+X(IndexSorted(C),:))/2));
                else
                    X(IndexSorted(i),:) = weight4(IndexSorted(i),:).*BestX(IndexSorted(i),:)-vb.*weight3(IndexSorted(i),:).*abs(tempPosition(A,:)-((X(IndexSorted(B),:)+X(IndexSorted(C),:))/2));
                end
            end
        end
        for i=1:(N/2)
            if  rand<VC2
                X(IndexSorted(i),:) = X(IndexSorted(i),j)*(1+randn(1));
            else
                A = randi([1,count]);
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    X(IndexSorted(i),:) = weight4(IndexSorted(i),:).*tempPosition(A,:)+vb.*weight3(IndexSorted(i),:).*abs(tempPosition(A,:)-((BestX(IndexSorted(i),:)+X(IndexSorted(i),:))/2)).*(1+Levy(dim));
                else
                    X(IndexSorted(i),:) = weight4(IndexSorted(i),:).*tempPosition(A,:)-vb.*weight3(IndexSorted(i),:).*abs(tempPosition(A,:)-((BestX(IndexSorted(i),:)+X(IndexSorted(i),:))/2)).*(1+Levy(dim));
                end
            end
        end 
    else
        shrink=2*(1-it/Max_iter); % a decreases linearly fron 2 to 0
        for i=1:size(X,1)
            if  rand<VC2
                X(i,:) = X(i,j)*(1+randn(1));
            else
                A = randi([1,count]);
                r = rand();
                vb = 2*shrink*r-shrink;%[-a,a]
                % Moving based on the bestPosition
                % The transformation range is controlled by weight3,bestPositions and X
                if r>VC1(i)
                    X(i,:) = weight4(i,:).*tempPosition(A,:)+vb.*weight3(i,:).*abs(0.5*(tempPosition(A,:)-BestX(i,:))).*(1+Levy(dim));
                else
                    X(i,:) = weight4(i,:).*tempPosition(A,:)-vb.*weight3(i,:).*abs(0.5*(tempPosition(A,:)-BestX(i,:))).*(1+Levy(dim));
                end
            end
        end  
    end
    
    
    Convergence_curve(it)=Destination_fitness;
    it=it+1;
end
toc
end



