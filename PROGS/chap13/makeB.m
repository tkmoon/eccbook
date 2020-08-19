
% Make the B matrix for an example code

outputs = [0 2 4 6
	       5 7 1 3
		   2 0 6 4
		   7 5 3 1];
dist = [0 0.6 2 3.4 4 3.4 2 0.6  % d(D0, x)
	    0.6 0 0.6 2 3.4 4 3.4 2  % d(D1, x)
	    2  .6 0 0.6 2 3.4 4 3.4   % d(D1, x)
	    3.4 2  .6 0 0.6 2 3.4 4   % d(D3, x)
	    4 3.4 2  .6 0 0.6 2 3.4   % d(D4, x)
	    3.4 4 3.4 2  .6 0 0.6 2   % d(D5, x)
	    2 3.4 4 3.4 2  .6 0 0.6   % d(D6, x)
	    0.6 2 3.4 4 3.4 2  .6 0];   % d(D6, x)

sp = '@{\kern 2pt}';
fout = fopen('bmat','w');
fprintf(fout,'\\left[\\begin{array}{l%sl%sl%sl%sl%sl%%\n',...
   sp,sp,sp,sp,sp);
fprintf(fout,'%sl%sl%sl%sl%%\n',sp,sp,sp,sp,sp);
fprintf(fout,'%sl%sl%sl%sl%sl%sc}\n',sp,sp,sp,sp,sp,sp);
for p=0:3
  for q=0:3
	% fprintf(fout,'%d%d: ',p,q);
	for p1=0:3
	  for q1=0:3
		out1 = outputs(p+1,p1+1);
		out2 = outputs(q+1,q1+1);
		d = dist(out1+1,out2+1);
		% print the stuff around the edges
		if(p1==0 & q1==0) % left edge stuff
		  fprintf(fout,'\\makebox[0cm][r]{\\makebox[3em][l]{%d%d}}',p,q);
		end
		if(p==0 & q==0) % top stuff
		  fprintf(fout,'\\rlap{\\smash{\\raisebox{1.5em}{\\makebox[.5em][c]{%d%d}}}}\n',p1,q1);
		end
		if(d==0)
		  fprintf(fout,'1');
		else
		  fprintf(fout,'x^{%g}',d);
		end
		if(~((p1==3) & (q1==3))) 
		  fprintf(fout,'&');
		else
		  if(~((p==3) & (q==3)))
			fprintf(fout,'\\\\');
		  end
		end
	  end
	end
	fprintf(fout,'\n');
  end
end
fprintf(fout,'\\end{array}\\right]');
fclose(fout);


% write out the p matrix
fout = fopen('pmat','w');
fprintf(fout,'\\left[\\begin{array}{l%sl%sl%sl%sl%sl%%\n',...
   sp,sp,sp,sp,sp);
fprintf(fout,'%sl%sl%sl%sl%%\n',sp,sp,sp,sp,sp);
fprintf(fout,'%sl%sl%sl%sl%sl%sc}\n',sp,sp,sp,sp,sp,sp);
firstrow = 1;
for p=0:3
  for q=0:3
	% fprintf(fout,'%d%d: ',p,q);
	if(p == q) continue; end;
	leftcol = 1;
	for p1=0:3
	  for q1=0:3
		if(p1==q1) 
		  if(p1==3) % last one
			fprintf(fout,'\\\\');
		  end
		  continue;
		end;
		out1 = outputs(p+1,p1+1);
		out2 = outputs(q+1,q1+1);
		d = dist(out1+1,out2+1);
		% print the stuff around the edges
		if(leftcol) % left edge stuff
		  fprintf(fout,'\\makebox[0cm][r]{\\makebox[3em][l]{%d%d}}',p,q);
		end
		if(firstrow) % top stuff
		  fprintf(fout,'\\rlap{\\smash{\\raisebox{1.5em}{\\makebox[.5em][c]{%d%d}}}}\n',p1,q1);
		end
		if(d==0)
		  fprintf(fout,'1');
		else
		  fprintf(fout,'x^{%g}',d);
		end
		if(~((p1==3) & (q1==3))) 
		  fprintf(fout,'&');
		else
		  if(~((p==3) & (q==3)))
			fprintf(fout,'\\\\');
		  end
		end
		leftcol = 0;
	  end
	end
	firstrow = 0;
	fprintf(fout,'\n');
  end
end
fprintf(fout,'\\end{array}\\right]');
fclose(fout);

% write out the m matrix
fout = fopen('mmat','w');
fprintf(fout,'\\left[\\begin{array}{l%sl%sl%sl%sl%sl%%\n',...
   sp,sp,sp,sp,sp);
fprintf(fout,'%sl%sl%sl%sl%%\n',sp,sp,sp,sp,sp);
fprintf(fout,'%sl%sl%sl%sl%sl%sc}\n',sp,sp,sp,sp,sp,sp);
firstrow = 1;
for p=0:3
  for q=0:3
	% fprintf(fout,'%d%d: ',p,q);
	if(p == q) continue; end;
	leftcol = 1;
	for p1=0:3
	  for q1=0:3
		if(p1~=q1) 
		  continue;
		end;
		out1 = outputs(p+1,p1+1);
		out2 = outputs(q+1,q1+1);
		d = dist(out1+1,out2+1);
		% print the stuff around the edges
		if(leftcol) % left edge stuff
		  fprintf(fout,'\\makebox[0cm][r]{\\makebox[3em][l]{%d%d}}',p,q);
		end
		if(firstrow) % top stuff
		  fprintf(fout,'\\rlap{\\smash{\\raisebox{1.5em}{\\makebox[.5em][c]{%d%d}}}}\n',p1,q1);
		end
		if(d==0)
		  fprintf(fout,'1');
		else
		  fprintf(fout,'x^{%g}',d);
		end
		if(~((p1==3) & (q1==3))) 
		  fprintf(fout,'&');
		else
		  if(~((p==3) & (q==3)))
			fprintf(fout,'\\\\');
		  end
		end
		leftcol = 0;
	  end
	end
	firstrow = 0;
	fprintf(fout,'\n');
  end
end
fprintf(fout,'\\end{array}\\right]');
fclose(fout);
