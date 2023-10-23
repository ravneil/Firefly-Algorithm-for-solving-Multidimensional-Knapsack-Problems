%
% 'Firefly Algorithm'
%

clc;
clear;
close all;

%% Problem Definition
for you=1:1%%%%User update
    global slooots;
    slooots= you;
 %storing in file
fileID = fopen('results_knapsack_FA_Modified_Final2023.txt','a+');
fprintf(fileID,'Problem %12.8f\r\n',slooots);

 Analysis_output=[];  
 Analysis_iteration=[]; 
  run=30;  
%Update me with number of runs needed.    
for me=1:run%%%%User update
% x = 0:.1:1;
% A = [x; exp(x)];
% c=[];
model = CreateModel(slooots);  % Create the Model

CostFunction = @(x) BinPackingCost(x, model);  % Objective Function

nVar = model.n%%%%2*model.n-1;     % Number of Decision Variables
VarSize = [1 nVar];     % Decision Variables Matrix Size

VarMin = 0;     % Lower Bound of Decision Variables
VarMax = 1;     % Upper Bound of Decision Variables


%% Firefly Algorithm Parameters

MaxIt=1000;         % Maximum Number of Iterations

nPop=50;            % Number of Fireflies (Swarm Size)
gamma=1;%1;%0.1;%            % Light Absorption Coefficient
beta0=2;%2;%1;%             % Attraction Coefficient Base Value
alpha=0.2;% 0.2;%1;%           % Mutation Coefficient%%%%%%%modified to 0.5 from 0.2
alpha_damp=0.98;    % Mutation Coefficient Damping Ratio

delta=0.05*(VarMax-VarMin);     % Uniform Mutation Range

m=2;


if isscalar(VarMin) && isscalar(VarMax)
    dmax = (VarMax-VarMin)*sqrt(nVar);
else
    dmax = norm(VarMax-VarMin);
end

%% Initialization

% Empty Firefly Structure
firefly.Position=[];
firefly.Cost=[];
firefly.Sol=[];

lim=0;
% Initialize Population Array
pop=repmat(firefly,nPop,1);
popC=repmat(firefly,nPop,1);
% Initialize Best Solution Ever Found
BestSol.Cost=0;%%%%inf;

% Create Initial Fireflies
for i=1:nPop
   pop(i).Position=round(randn(1,model.n));%%%%%%unifrnd(VarMin,VarMax,model.n);
   
    pop(i).Position=round(max(pop(i).Position,VarMin));
    pop(i).Position=round(min(pop(i).Position,VarMax));

   [pop(i).Cost, pop(i).Sol]=CostFunction(pop(i).Position);
  
   if pop(i).Cost>=BestSol.Cost%%%%%changed to > sicne max function problem <=BestSol.Cost
       BestSol=pop(i);
   end

end

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% Firefly Algorithm Main Loop

for it=1:MaxIt
    
    newpop=repmat(firefly,nPop,1);
    for i=1:nPop
        newpop(i).Cost = 0;%%%%%inf;
        newpop(i).Sol = 0;
          for j=1:nPop
            if pop(j).Cost >= pop(i).Cost %|| i==j

                 rij=norm(pop(i).Position-pop(j).Position)/dmax;
                beta=beta0*exp(-gamma*rij^m);
                e=delta*randn(VarSize); 

            if rij>0 || BestSol.Cost==1   
                
                newsol.Position = pop(i).Position ...
                                + beta*randn(VarSize).*(pop(j).Position-pop(i).Position);% ...
%                                  + rand;             
                newsol.Position=round(min(newsol.Position,VarMax));
                newsol.Position=round(max(newsol.Position,VarMin));

                [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position); %%%%%Evaluate
               
                if newsol.Cost >= newpop(i).Cost ||(newpop(i).Cost-newsol.Cost<10)
                    newpop(i) = newsol;
                    if newpop(i).Cost>=BestSol.Cost
                        BestSol=newpop(i);
                        
                    end
                end 
                
                if newsol.Cost < newpop(i).Cost && (newpop(i).Cost-newsol.Cost>10)
                pp=insertionFunction(newsol.Position,newsol.Cost, 1,model.n,model.limit); 
                newsol.Position=pp ;           
                   [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position); %%%%%Evaluate 
                   if newsol.Cost >= newpop(i).Cost %%%%%%%%%%%Update  
                    newpop(i) = newsol;
                    if newpop(i).Cost>=BestSol.Cost
                        BestSol=newpop(i);
                        
                    end
                   end
                end
            end
       if rij==0 || BestSol.Cost==1 

beta=rand; %1.5576;
                    newsol.Position = pop(j).Position ...%% change 2 - using pop Jth instead of ith term
                                + beta*randn(VarSize) ;%.*(pop(j).Position-pop(i).Position);%...
%                                      + alpha*e ; 

                    y = abs(tanh( 1*newsol.Position ));
                    s=rand;
                    y(y>=s)=1;
                    y(y<s)=0;    
                                
                    newsol.Position=y;
                    [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position); %%%%%Evaluate
   
                    if newsol.Cost >= newpop(i).Cost ||(newpop(i).Cost-newsol.Cost<5)
                        newpop(i) = newsol;
                        if newpop(i).Cost>=BestSol.Cost
                            BestSol=newpop(i);
                        end
                    end
                
                
                if newsol.Cost <= newpop(i).Cost &&(newpop(i).Cost-newsol.Cost>5)
                  pp=insertionFunction(newsol.Position,newsol.Cost,1,model.n,model.limit); 
                  newsol.Position=pp;               
                
                 [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position); %%%%%Evaluate 
                   if newsol.Cost >= newpop(i).Cost %%%%%%%%%%%Update  
                    newpop(i) = newsol;
                    if newpop(i).Cost>=BestSol.Cost
                        BestSol=newpop(i);
                        
                    end
                   end
                end
                end
%               end
            end
          end  
      
    end
    
    % Merge
    pop=[pop
         newpop];  %#ok    
    % Sort    
    [~, SortOrder]=sort([pop.Cost],'descend');
    pop=pop(SortOrder);    
    % Truncate
    pop=pop(1:nPop);
   
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    BestSol.Position;
    BestSol.Position(1);
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    temp_best_cost(it)=BestCost(it);   
    if it>10
        
%        steps=steps*alpha;
        if BestCost(it) == BestCost(it-9)
            
                break;
        end
    end
    % Damp Mutation Coefficient
    alpha = alpha*alpha_damp;

end

fprintf(fileID,'%6.2f %12.8f\r\n',me,BestCost(it));

Analysis_output(me)=BestCost(it);
Analysis_iteration(me)=it-9;


if (me==run)
 Mean_Output = mean(Analysis_output);
 Median_Output = median(Analysis_output);
 Max_Output = max(Analysis_output);
 Min_Output = min(Analysis_output);
 val = sum(Analysis_output == Max_Output);
 numc= sum(Analysis_output(:) == Max_Output);
 
 idx=find(Analysis_output==Max_Output)
 Best_Iteration = Analysis_iteration(idx(1));
 
  fprintf(fileID,'Min %12.8f\r\n',Min_Output);
 fprintf(fileID,'Median %12.8f\r\n',Median_Output);
 fprintf(fileID,'Max %12.8f\r\n',Max_Output);
 fprintf(fileID,'Mean %12.8f\r\n',Mean_Output);
 fprintf(fileID,'No. of Iteratation  %12.8f\r\n',Best_Iteration);
  fprintf(fileID,'Count Max %5.5f\r\n',numc); 

end

end

   
fclose(fileID);
BestSol.Position

%plot best convergence
figure;
plot(temp_best_cost,'LineWidth',2);
xlabel('Iteration');
ylabel('Fitness');

end

