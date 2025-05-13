function peaks = peakpick(s);
% peakpick: A peak picking algorithm
% 		Picks the peaks out in a vector by looking for
% 		negative zero crossings of the first derivative.
% 		peaks = peakpick(s)
%			peaks: a list of peaks found
%			s: the vector from which peaks are extracted.

% make 1st derivative
d1 = diff(s);

% find negative zero crossings (i.e. peaks)
[r,length] = size(d1);

peaks=[];

zerono = 1;
for l = 2:length
    if (d1(l-1)>0) & (d1(l)<=0)
       peaks(zerono) = l;
       zerono = zerono+1;   
    end;
end;
    


