function troughs = trofpick(s,t);
% trofpick: A trough picking algorithm
% 		Picks the troughs out in a vector by looking for
% 		positive zero crossings of the first derivative.
% 		troughs = trofpick(s)
%			troughs: a list of peaks found
%			s: the vector from which peaks are extracted.


if (nargin<2)
     tolerance = 0;
else
     tolerance = t;
end;

% make 1st derivative
d1 = diff(s);

% find negative zero crossings (i.e. peaks)
[r,length] = size(d1);

troughs=[];

zerono = 1;
for l = 2:length
    if (d1(l-1)<=tolerance) & (d1(l)>tolerance)
       troughs(zerono) = l;
       zerono = zerono+1;   
    end;
   end;

