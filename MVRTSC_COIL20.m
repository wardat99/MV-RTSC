clear;
clc
close all

addpath('dataset','functions');
load('COIL20_3VIEWS.mat')
A(:,:,1)=network_construction(X1',15);  % best at 15
A(:,:,2)=network_construction(X2',15);  % best at 15
A(:,:,3)=network_construction(X3',15);  % best at 15
num_views = 3;
truth=Y; clear Y X1 X2 X3
numClust = length(unique(truth));


mu = 9;              %best at  9  
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