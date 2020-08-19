
% compute probabilities for the quantization of the 
% Euclidean metric (for Viterbi decoding in Gaussian noise)

sigma = 1;

a1 = 1;
a2 = -1;								% signal amplitudes

q1 = -1;								% quantization levels
q2 = 0;
q3 = 1;

% compute the probability of falling in the quantized region
p4 = qf((q3-a1)/sigma);
p3 = qf((q2-a1)/sigma) - p4;
p2 = qf((q1-a1)/sigma) - (p3+p4);
p1 = 1-(p2+p3+p4);

% compute log probabilities
p4l = -log(p4);
p3l = -log(p3);
p2l = -log(p2);
p1l = -log(p1);


% subtract off the smallest probability
p4d = 0;
p3d = p3l-p4l;
p2d = p2l-p4l;
p1d = p1l-p4l;


% search for a scale factor that will result in smallest truncation error
minm = 2;
clist = 1:.001:5;
maxpd = max([p1d,p2d,p3d,p4d]);
maxqlist = [];
mlist = [];
clist2 = [];
maxdespd = 8;
for c = clist
  r4 = p4d*c - floor(p4d*c);
  r3 = p3d*c - floor(p3d*c);
  r2 = p2d*c - floor(p2d*c);
  r1 = p1d*c - floor(p1d*c);
  if(c*maxpd > maxdespd)  % if we have the largest value we want to get
	break;
  end
  maxqlist = [maxqlist c*maxpd];
  clist2 = [clist2 c];
%  m = max([r4 r3 r2 r1]);
  m = sum([r4 r3 r2 r1]);
  mlist = [mlist m];
  if(m < minm) 
	cmin = c;
	minm = m;
  end
end
plot(clist2,mlist); hold on;
plot(clist2,maxqlist);

% the un-truncated probabilities
cmin
p4d*cmin
p3d*cmin
p2d*cmin
p1d*cmin
