% Color pallette for different types for making figures with.

rationalisedTypeNames = {'','PL','PLN','ChS','ChT','On','PBU','Unc'};
rationalisedTypeNums = [0:7];

typecolormap = [ ...
                arrayfun(@(x) hex2dec(x),'05e')/16; ...         % PL
                arrayfun(@(x) hex2dec(x),'2af')/16; ...         % PLN
                arrayfun(@(x) hex2dec(x),'f33')/16; ...         % ChS
                arrayfun(@(x) hex2dec(x),'fa0')/16; ...         % ChT
                arrayfun(@(x) hex2dec(x),'3a2')/16; ...         % On
                arrayfun(@(x) hex2dec(x),'859')/16; ...         % PBU    
                arrayfun(@(x) hex2dec(x),'045')/16; ...         % Unc
                ];

            % Function to find the "colour number"          
findRatTypeNum = @(x)  rationalisedTypeNums( find( strcmp( x,  rationalisedTypeNames ) ) );    
findRatTypeCol=  @(x) typecolormap(findRatTypeNum(x),:);
