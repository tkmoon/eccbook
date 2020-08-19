function writesparse(A,fp)
% function writesparse(A,fp)
%
% Write the binary matrix A into the file fp
% in MacKay's sparse format

% Todd Moon

[m,n] = size(A);
fprintf(fp,'%d %d\n',n,m);   % dimensions of A
fprintf(fp,'%d %d\n',max(sum(A,1)),max(sum(A,2)));  % max column and row weights
fprintf(fp,'%d ',sum(A,1)); fprintf(fp,'\n');  % column weights (array)
fprintf(fp,'%d ',sum(A,2)); fprintf(fp,'\n');  % row weights (array)
for i=1:n  % write out the nonzero positions in each column
  fprintf(fp,'%d ',find(A(:,i)));
  fprintf(fp,'\n');
end
for i=1:m  % write out the nonzero positions in each row
  fprintf(fp,'%d ',find(A(i,:)));
  fprintf(fp,'\n');
end