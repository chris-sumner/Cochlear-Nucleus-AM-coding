function output = NHPPsurrogate(spikePeriods, per, ISI, NBrepe, refrac, win,PHrate,forceZeroRefrac)
% Return a surrogate spike train using
%
% spikePeriods is the period assigned to each spike time
% per is period in ms or 1000/fm

if nargin<8
    forceZeroRefrac = 0;
end

if per<0.2 % I doubt fm will be large enough for 1000/fm to satisfy this
    error('This will go a lot faster if you express the period in ms');
end

winDur = win(2)-win(1); % duration of surrogate spike train

NBspk = length(spikePeriods); % 1 period for each spike so this gives total number of spikes
ISIcell_PHsur = cell(1,NBrepe);

% rate = spikes per ms
spikeRatePerMs = NBspk/(NBrepe*winDur);       % Window duration (100 ms ) -> spikesPerSweep in kHz

% if the refractory period is greater than half over the number of spikes per ms then
% limit the refractory period to half of the spike rate per ms - WHY?
if refrac>(0.5/spikeRatePerMs)
    refrac = 0.5/spikeRatePerMs;
    fprintf('Warning: refactory period limited to half of the spikeRatePerMs\n');
end;

% this bit checks that the refractory time is big enough - WHY?
Rd = spikeRatePerMs/(1-spikeRatePerMs*refrac); % maybe this should be refrac - spikeRate * refrac?????
if(Rd<0)
    disp('Refractory time too small - limited to minimum ISI')
    % limit refrac to smallest ISI
    refrac = min(ISI);
    Rd = spikeRatePerMs/(1-spikeRatePerMs*refrac);
    if(Rd<=0)
        disp('Error in NHPPsurrogate at line 30: Rd is negative!')
        st_sur=[]; ISI_sur=[]; ph_sur=[]; ISIcell_PHsur=[];
        output.st = st_sur;
        output.ISI = ISI_sur;
        output.ph = ph_sur;
        output.ISIcell = ISIcell_PHsur;
        output.refrac = refrac;
        return;
    end
end

if forceZeroRefrac
    refrac = 0;
end

% Get a waveform of the NHPP process intensity
% First create a ph of the spikePeriods and normalise by spikes per bin
dPH = 0.02;
dPH2 = 0.002;
% if PHrate is not submitted as an input argument
if nargin<7
    PHrate = histc(spikePeriods,0:dPH:1);
end
PHrate = PHrate/mean(PHrate);
PHrate(end)=PHrate(1); % to wrap around and interpolate?
% this bit seems to create an interpolated copy of the ph with a much finer
% resolution (dPH2 is 10 times dPH), normalised to the Rd
sPHrate = Rd*interp1(0:dPH:1,PHrate,0:dPH2:1-dPH2,'cubic');

% ***** This carries an extra loop to estimate the error in rate *****
% I believe we have the extra loop so that we can increase the length of
% the simulations duration and thus estimate rate with more observations
dt = dPH2*per; % finer rate x period
modRate = repmat(sPHrate,1,floor(winDur/per)); % for full window duration
cpt = 1;
beg = win(1);
simDur = 5*winDur; % simulation duration -> longer = more reliable measure of rate
while(beg<(win(1)+simDur))
    rateLth = length(modRate);
    t = beg:dt:(rateLth*dt+beg-dt);
    U = rand;
    E = -log(U);
    idx2shift = floor(beg/dt);
    
    % CIRCSHIFT REPLACED - SEE whyDoesNHPPtakeSoLong.m
    %shiftRate = circshift(modRate,[0 -idx2shift]);
    newshift = rem(idx2shift,length(modRate));
    shiftRate=[modRate(newshift+1:end) modRate(1:newshift)];
    
    y = cumtrapz(heav(t-beg-refrac).*shiftRate)*dt; 
        % differential equation integration - RIP Jonathan. You were amazing. 
    maxy = max(y);
    
    % this code replaces the while loop and should ensure that we
    % select a random number that results in maxy being greater
    % than E
    maxRandom = exp(-maxy);
    U=maxRandom+rand*(1-maxRandom); % 1-maxRandom gives the range of U that will give value E less than ymax
    E = -log(U);
    newtime = find(y>=E,1,'first')*dt;
    cpt = cpt+1;
    st4rate(cpt) = beg + newtime;
    beg = st4rate(cpt);
end

simRate = length(st4rate)/simDur;       % rate in kHz

dt = dPH2*per;
newError = min(simRate, 0.99); % why did I put this in - because errors need to be less than 1 to enter the following loop
newError = max(newError, 0.051); % this ensures that the following while loop is always entered on the first iteration
c = 1;
% error = diff(input rate, simulated rate) / input rate  
while (abs(newError)>0.05 && c<12 && abs(newError)<1) % while error is between 0.05 and 1 and repeats < 12
    sPHrate = sPHrate*(1+newError); % new error has to be between 0 and 1 for this to not be crazy methinks
    modRate = repmat(sPHrate,1,floor(winDur/per));
    prevERR = newError;
    st_cell = {};
    st_sur = [];
    isi_sur = [];
    %     isiCnt = 1;
    for ii=1:NBrepe
        cpt = 1;
        beg = win(1);
        while(beg<win(2)) % goes through window from start to finish
            rateLth = length(modRate);
            t = beg:dt:(rateLth*dt+beg-dt);
            idx2shift = floor(beg/dt);
            
            % CIRCSHIFT REPLACED - SEE whyDoesNHPPtakeSoLong.m
            %shiftRate = circshift(modRate,[0 -idx2shift]);
            newshift = rem(idx2shift,length(modRate));
            shiftRate=[modRate(newshift+1:end) modRate(1:newshift)];
            
            y = cumtrapz(heav(t-beg-refrac).*shiftRate)*dt;
            %             counter
            maxy = max(y);
            
            % this code replaces the while loop and should ensure that we
            % select a random number that results in maxy being greater
            % than E
            maxRandom = exp(-maxy);
            U=maxRandom+rand*(1-maxRandom); % 1-maxRandom gives the range of U that will give value E less than ymax
            E = -log(U);
            newtime = find(y>=E,1,'first')*dt;
            st_cell{ii}(cpt) = beg + newtime;
            beg = st_cell{ii}(cpt);
            cpt = cpt+1;
        end
        st_temp = [st_cell{ii}];
        st_temp = st_temp( st_temp>win(1) & st_temp<win(2) ); % only keep spiketimes in the analysis window
        ISIcell_PHsur{ii} = diff(st_temp);
        st_sur = [st_sur st_temp];
        isi_sur = [isi_sur diff(st_temp)];
    end
    st_sur = st_sur(:);
    % error term is used to try to force the model output to have a similar
    % spikerate per ms to that specified by the input variables
    newRate = length(st_sur)/(NBrepe*winDur); % in spikes per ms
    newError = (spikeRatePerMs-newRate)/spikeRatePerMs;
    disp(['Error after spikeRatePerMs (' num2str(c)  '):' num2str(newError)])
    % we only ever take the spiketrain with the lowest error
    if abs(newError)<abs(prevERR) || c==1
        best_st_sur = st_sur;
        best_isi_sur = isi_sur;
        best_ISIcell_PHsur = ISIcell_PHsur;
    end
    c = c+1;
end

st_sur = best_st_sur;
ISIcell_PHsur= best_ISIcell_PHsur;
ISI_sur = best_isi_sur;

ip_sur = floor(st_sur/per);
ph_sur  = (st_sur - ip_sur*per)/per;

output.st = st_sur;
output.ISI = ISI_sur;
output.ph = ph_sur;
output.ISIcell = ISIcell_PHsur;
output.refrac = refrac;



function y = heav(x)
y = 1.0*(x>0);