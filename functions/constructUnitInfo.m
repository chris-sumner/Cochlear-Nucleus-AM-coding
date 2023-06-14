function  unitinfo = constructUnitInfo(wholeList171212,EXPLOGLIST,typeNames,thisdatastartind,thisdatainds,newdatainds,conditionswithspikes)
% Make unitinfo structure. 

% Get the am frequencies, VS, Z scores
allamfreqs = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'modFreq')));
VSlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'VS1')));  %30);
Rayleighlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'Ray1')));  %32);
Zlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'newZs')));  %66);
MLtest = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'newZaboveRefrac1')));  %70);
SperP = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'spikesPerPeriod')));  %17);
ref0s = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'newZaboveRefrac0')));  %69);
cilist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'CI')));  %60);
ciPlist = wholeList171212(thisdatainds,find(strcmp(EXPLOGLIST,'CIpVal')));  %71);

% A load of useful information about the data is stored in unitinfo.
unitinfo.unitid = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'uniqueID')));   % dividingdata(newdatainds(thisdatastartind),1);
unitinfo.modLevel = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'modLevel')));   % dividingdata(newdatainds(thisdatastartind),2);
unitinfo.amToneDur = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'amToneDur')));   % dividingdata(newdatainds(thisdatastartind),3);
unitinfo.carrFreq =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'carrFreq')));   % dividingdata(newdatainds(thisdatastartind),4);
unitinfo.depthMod =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'depthMod')));   % dividingdata(newdatainds(thisdatastartind),5);
unitinfo.cf =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'cf')));   
unitinfo.cf_thr =  wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'cf_thr')));   
unitinfo.allamfreqs =  allamfreqs;
unitinfo.usedamfreqs =  allamfreqs(conditionswithspikes);
unitinfo.TypeNum = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'TypeNum')));  %3);
unitinfo.TypeName = typeNames{wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'TypeNum')))+1};
unitinfo.VS = VSlist;
unitinfo.Rayleigh = Rayleighlist;
unitinfo.Z = Zlist;
unitinfo.ref0s = ref0s;
unitinfo.ci = cilist;
unitinfo.cipVal = ciPlist;
unitinfo.MLPoissTest = MLtest;
unitinfo.CV = wholeList171212(newdatainds(thisdatastartind),find(strcmp(EXPLOGLIST,'newCV')));  %63);
unitinfo.spikesperperiod = SperP;
unitinfo.conditionswithspikes = conditionswithspikes;
unitinfo.wholeListInds = thisdatainds;

% This keeps the entire set of statistics for the FIRST condition.
unitinfo.scholesWholeListFirstEntry = wholeList171212(newdatainds(thisdatastartind),:);
    