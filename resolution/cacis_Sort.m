function z = cacis_Sort(Z,positions,m,n,p)
% Z vec(Rxx)
%positions: 阵元坐标
%m
%n
%p：压缩系数；
s = 1;
row_num = zeros(2*m*n-2*m*(n-1)/p-1,1);
row_matrix = zeros(m+n-1,m+n-1);
for i = 1:m+n-1
    row_matrix(i,:) = positions(i) - positions;
end
for j = -m*n+m*(n-1)/p+1:m*n-m*(n-1)/p-1
    a = find(row_matrix == j);
    row_num(s) = a(1);
    s = s+1;
end
    row_num = row_num';
    z = Z(row_num);
end
        
