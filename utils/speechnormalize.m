function Y = speechnormalize(X)
[m,n] = size(X);
for i =1:n
     A(1,i)=max(X(:,i));
%    A(1,i)=norm(X(:,i));
end
A = repmat(A, m, 1);
Y=X./A;
end