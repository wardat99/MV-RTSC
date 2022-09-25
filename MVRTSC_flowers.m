clear
close all
clc

addpath('dataset','functions');

load('flowers.mat');
num_views = length(data);
numClust = length(unique(truth));
%% 

sigma=5; % best at 5
A(:,:,1)=data{1}./(2*sigma^2); 
A(:,:,2)=data{2}./(2*sigma^2); 
A(:,:,3)=data{2}./(2*sigma^2); 
[n1,n2,n3]=size(A);
for v=1:1:num_views
    for i=1:n1
     for  j=1:n2
     A(i,j,v)=exp(-1*A(i,j,v));  % create adjacency tensor 
     end
    end 
end 

mu = 4;              %best at  4  
lambda =  0.1;       %best at  0.1


for i=1:1
    fprintf('----MV-RTSC start, attempt number %d--------\n', i);
[Plabel,Timecost(i)] = MVRTSC(A,mu,lambda,numClust);
   fprintf('----MV-RTSC end, attempt number %d--------\n', i);
acc(i) =  Compute_accuracy(truth,Plabel);
        [Aa nmi(i) avgenti] = compute_nmi(truth,Plabel);
        [f(i),p(i),r(i)] = compute_f(truth,Plabel);
        if (min(truth)==0)
            [AR(i),RI,MI,HI]=RandIndex(truth+1,Plabel);
        else
            [AR(i),RI,MI,HI]=RandIndex(truth,Plabel);
        end 
end
fprintf('Acc: %.2f  (%.2f)\n' , mean(acc), std(acc));
fprintf('nmi: %.4f  (%.4f)\n' , mean(nmi), std(nmi));
fprintf('AR: %.4f   (%.4f)\n' , mean(AR), std(AR));
fprintf('F: %.4f    (%.4f)\n' , mean(f), std(f));
fprintf('P: %.4f    (%.4f)\n' , mean(p), std(p));
fprintf('R: %.4f    (%.4f)\n' , mean(r), std(r));
fprintf('Timecost: %.4f  \n\n' , mean(Timecost));