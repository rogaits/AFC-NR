function [W,X] = mwf_GEVD2(Ryy,Rnn,Rank)
P = Rank-1;
[N,~] = size(Ryy);
[X,Lq]=eig(Ryy,Rnn);
[sort_values, sort_index] = sort(diag(Lq),'descend');
Lq = diag(sort_values);            
X = X(:,sort_index); 
Q = inv(X');
if Rank==1
    Q = Q(:,1);
    D = 1 - 1./diag(Lq); % Rank-1 approximation based on eigenvalues    
    D = diag(D(1));
    W = X(:,1)*D*Q';
%     W = W(:,[1 2]); % first column is for the loudspeaker signal
else
    Q = Q(:,1:Rank);
    D = 1 - 1./diag(Lq); % Rank-1 approximation based on eigenvalues    
    D = diag(D(1:Rank));
    W = X(:,1:Rank)*D*Q';
%     W = W(:,[1 2]); % first column is for the loudspeaker signal
end
% Sigma_ss = X'*(Ryy-Rnn)*X;

%  P+2 is the first microphone signal
end