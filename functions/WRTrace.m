function op = WRTrace(S,T,a,Tlim)

% Last spike. 
maxT_ms = max(Tlim);
minT_ms = min(Tlim);

S = S(S>=minT_ms &  S<maxT_ms);
T = T(T>=minT_ms &  T<maxT_ms);

% Make time steps
dt   = 0.05;
allT = maxT_ms;
ts   = 0:dt:allT;
N    = length(ts);

% Reserve space for spike train
yS = zeros(N, 1);
yT = zeros(N, 1);

% Make spike train
yS(round(S / dt)+1) = 1;
yT(round(T / dt)+1) = 1;

op.S = yS;
op.T = yT;

% Filter waveforms
as = [1 -1+a*dt];
as = conv(as, as);
yS = filter(dt, as, yS);
yT = filter(dt, as, yT);

op.yS = yS;
op.yT = yT;
op.Tlim = Tlim;

op.diffSqTrace = (yS - yT).^2;
op.sumDiffSq = sum( op.diffSqTrace );
