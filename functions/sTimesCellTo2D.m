function spikesOut = sTimesCellTo2D(spikeCell,nSweeps,opt)

% converts a cell array / 1D array of spike times in to a 2D array with sweep x spike
% times, NaNs for all the empty values

if nargin<3
    opt=1;
else
    opt=2;
end

switch opt
    case 1
        % make sure input is the right dimensionality
        spikeCell = spikeCell(:)';
        
        if iscell(spikeCell)
            cellLengths = cellfun('length',spikeCell);
            spikesOut = NaN*ones(max(cellLengths),nSweeps);
            for ii = 1:nSweeps
                spikesOut(1:cellLengths(ii),ii) = spikeCell{ii}';
            end
        else
            sweepInds = [0 find(diff(spikeCell)<0) length(spikeCell)];
            sweepDiffs = diff(sweepInds);
            spikesOut = NaN*ones(max(sweepDiffs),nSweeps);
            for ii = 1:length(sweepDiffs)
                spikesOut(1:sweepDiffs(ii),ii) = spikeCell(sweepInds(ii)+1:sweepInds(ii+1))';
            end
        end
        % 2D array to cell
    case 2
        for ii = 1:nSweeps
            theseVals = spikeCell(:,ii);
            spikesOut{ii} = theseVals(~isnan(theseVals))';
        end
        spikesOut = spikesOut(:);
end