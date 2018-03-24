function z = CADiS_Sort(Z,positions,m,n,p)
s = 1;
row_matrix = zeros(m+n-1,m+n-1);

for i = 1:m+n-1
    row_matrix(i,:) = positions(i) - positions;
end
if m/p == 1
    aa = m*n+m/p-1;
    bb = -m*n-m/p+1;
    row_num = zeros(1,aa-bb+1);
    for i = bb:aa
        a = find(row_matrix == i);
        row_num(s) = a(1);
        s = s+1;
    end
else
    aa = m*n+m/p-1;
    bb = (m/p-1)*(n-1);
    row_num = zeros(1,aa-bb+1);
    for k = bb:aa
        a = find(row_matrix == k);
        row_num(s) = a(1);
        s = s+1;
    end
end
    row_num = row_num';
    z = Z(row_num);
end
