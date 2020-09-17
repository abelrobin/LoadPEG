function [BPC] = basePeakChrom(specdata,masses,eicGrid)

%Please attribute Robin J. Abel in resulting works.
%
%Current Version 1.0, 2019-11-01
%
%A function to generate a base peak chromatogram from full scan specdata
%
%Usage - [BPC] = basePeakChrom(specdata)
%
%   specdata	full scan data from a chromatogram
%   eicGrid		matrix of row vectors to form EICs -> BPC if desired

%% Version History
% Version 1.0, 2019-11-01, first version

%Generate base peak chromatogram
if exist('eicGrid','var')
    specdata2 = zeros(size(specdata,1),size(eicGrid,1));
    for i = 1:size(eicGrid,1)
        currEic = eicGrid(i,:);
        currEic = currEic(currEic ~= 0);
        [eicLocs,~] = ismember(masses,currEic);
        specdata2(:,i) = sum(specdata(:,eicLocs),2);
    end
    specdata = specdata2;
end
BPC = max(specdata,[],2);
end