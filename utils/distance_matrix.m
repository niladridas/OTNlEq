% This function gives the distance values in form of a matrix D
function D = distance_matrix(X_f)
M = size(X_f,2);
D = zeros(M,M);
for i = 1:M-1
    for j = i+1:M
        D1 = norm(X_f(:,i)-X_f(:,j));
        D(i,j) = D1;
        D(j,i) = D(i,j);
    end
end
end