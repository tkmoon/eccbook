function binv = invmodp(b, p)
% function binv = invmodp(b, p)

if(b == 0)
  error('Error: Illegal division by zero');
end
[g,x,y] = gcd(p,b);
if(g ~= 1)
  error('Error: attempting to divide by noninvertible number');
end
binv =  mod(y,p);

