%
% Fitness function
%

function [z, sol] = BinPackingCost(q, model)

    n = model.n;%%number of variables
    c = model.c;%number of knapsack
    price_set=model.price_set;
    weight_set=model.weight_set;%%m resources with capacities
    limit=model.limit;%%knapsack limit

sum_of_prices=[];
sum_of_weights=[];
q;



for i= 1:c
       
        sum_of_prices(i)=sum(q.*price_set(i,:));
        sum_of_weights= sum(q.*weight_set);          %q.* weight_set[i][j]);
       
end   
sum_of_prices;
limit;
        if sum_of_prices <= limit   %625
            z= sum_of_weights;
        else
            
            z = 1;% penalty function
        end

    sol.nBin = weight_set;
    sol.B = price_set;
   
    
end