% Set up parameters for the irregular repeat/accumulate codes used
% in DVB standards

% rate = n/m;
nrate = 2;
mrate = 3;
[ira,q,kldpc] = getiradat(nrate,mrate);
Nldpc = 64800;

Nminusk = Nldpc - kldpc;

nblockperrow = 360;    %from the standard, 302.755, p. 37

H1 = zeros(Nminusk,Nldpc,'uint8');

for i1 = 1:length(ira)  % for each row of ira data
   baseinfoidx = (i1-1)*nblockperrow;
   row = ira{i1};
   for m=0:nblockperrow - 1
      for i2 = 1:length(row)
         x = row(i2);
         % fprintf('i1=%d  x=%d ',i1,x);
         infoidx = baseinfoidx + m;
         paridx = mod(x + mod(m,360)*q,Nldpc - kldpc);
         % fprintf('p%d = i%d (+) %d\n',paridx,infoidx,paridx);
         if(H1(paridx +1,infoidx +1))  % already set
            fprintf('Bit already set at H(%d,%d)\n',paridx, infoidx);
         end            
         H1(paridx +1,infoidx +1) = H1(paridx +1,infoidx +1) + 1;
         % fprintf('H(%d,%d)+=1\n',paridx,infoidx);
      end
   end
end

% Add in the columns of H1 coresponding to the parity bits
for i=1:Nminusk
   H1(i, kldpc+i) = 1;
   % fprintf('H(%d,%d)=1\n',i-1,kldpc+i-1);
end

% Add in the columns of H1 corresponding to the final accumulation
for i=1:Nminusk - 1
   % fprintf('p%d = p%d (+) p%d\n',i,i,i-1);
   H1(i +1,kldpc + i) = 1;
end

% extract the sparse representation

% row information
maxrowwt = 0;
for m = 1:Nminusk
   N{m} = find(H1(m,:)) - 1;    % zero-based indexing
   maxrowwt = max([maxrowwt length(N{m})]);
end

% column information
maxcolwt = 0;
for n=1:Nldpc
   M{n} = find(H1(:,n)) - 1;    % zerobased indexing
   maxcolwt = max([maxcolwt length(M{n})]);
end

fp = fopen(sprintf('DVBT2LDPCR%d%d.txt',nrate,mrate),'w');
fprintf(fp,'%d %d\n',Nldpc,Nminusk);
fprintf(fp,'%d  %d\n',maxcolwt,maxrowwt);
for n=1:Nldpc
   fprintf(fp,'%d ',length(M{n}));
end
fprintf(fp,'\n');
for m=1:Nminusk
   fprintf(fp,'%d ',length(N{m}));
end
fprintf(fp,'\n');

for n=1:Nldpc
   fprintf(fp,'%d ',M{n});
   fprintf(fp,'\n');
end
for m=1:Nminusk
   fprintf(fp,'%d ',N{m});
   fprintf(fp,'\n');
end
fclose(fp);



H2 = double(H1);


