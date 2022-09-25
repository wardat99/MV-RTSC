function Gs = getGuassianSimilarity(x,y,sigma)
 
 %---------------
 
 % Mohammad Alwardat
 
 %-----------------
if isvector(x)==0 || isvector(y)==0
    error('x and y have to be vectors!')
end
if length(x)~=length(y)
    error('x and y have to be same length!')
end
xy   = x-y;
sigma_sq = (sigma)^2;
n_sq   = (norm(xy))^2;
z = (n_sq)/(2*sigma_sq);
 if z == inf
     z = 0;
 end
Gs   = exp (-z);
