function [A1] = network_construction(X,sigma)


    [n,m] = size(X);
   
    for i=1:m
        for j=1:m
            A1(i,j) = getGuassianSimilarity(X(:,i),X(:,j),sigma);
        end 
    end


end

