function [pegstruct] = LoadPEG(pegFile, procspec)
%Please read LICENSE.md before using, and attribute Dr. Robin J. Abel and Dr. James J. Harynuk
%with appropriate citation using DOI: 10.5281/zenodo.4035154
%
%Current Version 1.9, 2020-09-17
%
%LoadPEG.m		Loads LECO .peg files (up to the final version of ChromaTOF v4.x) into a Matlab
%				struct. Requires the PEGUnpack.mexw64 file for decoding of spectral data. Mac
%				version available upon request.
%
%Usage: [pegstruct] = LoadPEG(pegFile, procspec)
%
% inputs:
%	 pegFile				Path to .peg file as a string or character array
%	procspec				If 0, don't process chromatographic data
%
% returns:
%	pegstruct.fileName		Name and location of imported file
%	pegstruct.time			Centred timestamps vector, in milliseconds
%	pegstruct.time1D		1D timestamps in seconds
%	pegstruct.numScans		Number of scans
%	pegstruct.startMass		Lowest mass
%	pegstruct.massRange		Number of masses
%	pegstruct.masses		Index of specdata column m/z values
%	pegstruct.dataRate		Acquisition data rate
%	pegstruct.acqDate		Acquisition date string
%	pegstruct.modTime		Modulation period, in seconds (0 for 1D)
%	pegstruct.acqDelay		Delay before first scan, in seconds (0 for no delay)
%	pegstruct.specdata		Spectrum matrix[scan, mass]
%	pegstruct.tic			Total mass intensity vector by scan
%	pegstruct.bpc			Base peak intensity vector by scan

%% Version History
% Version 1.0, around 2018-07-01, (history lost)
% Version 1.1, 2018-08-27
% Version 1.2, around 2018-09-01
%	added offsets for loading in modulation period and solvent delay
% Version 1.3, around 2018-11-20
%	wrote .mex code to fast unpack spectra using OpenMP parallelization
% Version 1.4, around 2019-01-05
%	corrected time vector due to 0-indexing problem after code translation
% Version 1.5, 2019-06-11
%	added support for importing variable modulation files and recording of
%	the version of ChromaTOF used for acquisition
% Version 1.6, 2019-06-27
%	refined the modulation period seek position from fFR@ to a preceding
%	hex sequence which appears to be more stable
% Version 1.7, 2019-11-16
%	added pegstruct.bpc to output structure
% Version 1.8, 2020-01-27
%	added offsets to gather MS parameters for use when converting to netCDF format
% Version 1.9, 2020-09-17
%	first version for publication on GitHub

% Zero scalar return values
startMass = 0;
massRange = 0;
dataRate = 0;
numScans = 0;

if ~exist('procspec', 'var')
	procspec = 1;
end

% Open the file
if ~exist(pegFile, 'file')
	errmsg = strcat(pegFile, ' not found, check input path and name');
	error(errmsg)
else
	fid = fopen(pegFile);
end

% Check file open
if fid > 2
	% Read the magic number equivalent to the strings Hs&S for Peg4 or Gs&S for earlier formats
	magic = fread(fid, 1, '*int');
	if (magic ==  1395028808 || magic ==  1395028807)
		% Read the section offsets
		offsets = fread(fid, 6, '*int');

		% Read the last double in the file and convert to Matlab formatted date
		fseek(fid, offsets(6)-8, 'bof');
		datedays = 693960 + fread(fid, 1, 'double');
		acqDate = datestr(datedays);

		% Read the header as a byte block
		fseek(fid, offsets(2), 'bof');
		headerData = fread(fid, offsets(3)-offsets(2), '*uint8');

		% Find the version of ChromaTOF used for acquisition
		currSegmentOffset = 0;
		segmentFlag = 'Version';
		flagData = bin2dec(dec2bin(segmentFlag));
		for i = 1:(offsets(3)-(offsets(2)+11))
			% Section starts with the string Version
			if headerData(i:i+(length(segmentFlag)-1)) ==  flagData
				% Store the offset and stop search
				currSegmentOffset = offsets(2) + i + length(segmentFlag)-1;
				break;
			end
		end
		if (currSegmentOffset > 0)
			% Read & overshoot the version text, extract number using regular expression
			fseek(fid, currSegmentOffset-8, 'bof');
			acqVer = char(regexp(fread(fid, 64, '*char')', ' ([0-9.])+ ', 'tokens', 'once'));
		end

		% Find the Field Service section of the header
		currSegmentOffset = 0;
		segmentFlag = 'FIELD SERVICE';
		flagData = bin2dec(dec2bin(segmentFlag));
		for i = 1:(offsets(3)-(offsets(2)+11))
			% Section starts with the string CRampStepd
			if headerData(i:i+(length(segmentFlag)-1)) ==  flagData
				% Store the offset and stop search
				currSegmentOffset = offsets(2) + i + length(segmentFlag)-1;
				break;
			end
		end

		% Pull mass spec parameters from header
		fseek(fid, currSegmentOffset+331, 'bof');
		sourceTemp = double(fread(fid, 1, 'double'));
		fseek(fid, 196, 'cof');
		ionizVolt = double(fread(fid, 1, 'double'));
		fseek(fid, 80, 'cof');
		accelVolt = double(fread(fid, 1, 'double'));
		fseek(fid, 64, 'cof');
		tuneVolt = double(fread(fid, 1, 'double'));
		fseek(fid, 8, 'cof');
		transLineTemp = double(fread(fid, 1, 'double'));
		msParams.sourceTemperature = sourceTemp;
		msParams.ionizationVoltage = ionizVolt;
		msParams.accelerationVoltage = accelVolt;
		msParams.tuneVoltage = tuneVolt;
		msParams.transferLineTemperature = transLineTemp;
		%%%CAUTION: the offsets for the above values have not been fully confirmed to be correct%%%

		% Find the CItrSegment section of the header
		currSegmentOffset = 0;
		segmentFlag = 'CItrSegment';
		flagData = bin2dec(dec2bin(segmentFlag));
		for i = 1:(offsets(3)-(offsets(2)+11))
			% Section starts with the string CItrSegment
			if headerData(i:i+(length(segmentFlag)-1)) ==  flagData
				% Store the offset and stop search
				currSegmentOffset = offsets(2) + i + length(segmentFlag)-1;
				break;
			end
		end

		if (currSegmentOffset > 0)
			% Read from the CItrSegment header segment
			fseek(fid, currSegmentOffset + 18, 'bof');
			startMass = double(fread(fid, 1, '*uint16'));
			massRange = double(fread(fid, 1, '*uint16'));
			dataRate = double(fread(fid, 1, 'double'));

			% Read from the spectrum statistics block
			fseek(fid, offsets(4) + 8198, 'bof');
			numScans = double(fread(fid, 1, '*int'));

			if procspec ==  1
				% Read the spectral data as unsigned short
				fseek(fid, offsets(3) + 4, 'bof');

				% Call .mex to perform rapid bit unpacking
				specdata = PEGUnpack(fread(fid, [numScans, massRange], '*uint16'));

				% Calculate TIC (sum all masses for each scan)
				tic = sum(specdata, 2);
			end

			% Calculate scan-centered time vector in ms
			time(1:numScans) = 1000/dataRate;
			toffset = 500/dataRate;
			for i = 1:numScans
				time(i) = toffset + (time(i)*(i-1));
			end
		end

		% Find the CRampStep section of the header
		currSegmentOffset = 0;
		segmentFlag = 'CRampStepd';
		flagData = bin2dec(dec2bin(segmentFlag));
		for i = 1:(offsets(3)-(offsets(2)+11))
			% Section starts with the string CRampStepd
			if headerData(i:i+(length(segmentFlag)-1)) ==  flagData
				% Store the offset and stop search
				currSegmentOffset = offsets(2) + i + length(segmentFlag)-1;
				break;
			end
		end

		% Find the value 0000000001006400 after CRampStepd
		subSegmentOffset = 0;
		segmentFlag = ['00';'00';'00';'00';'01';'00';'64';'00'];
		flagData = hex2dec(segmentFlag);
		for i = 1:(offsets(3)-currSegmentOffset)
			% Section starts with the hex series 0000000001006400 following CRampStepd
			if headerData(i:i+(length(segmentFlag)-1)) ==  flagData
				% Store the offset and stop search
				subSegmentOffset = offsets(2) + i + length(segmentFlag)-1;
				break;
			end
		end
		if (subSegmentOffset > 0)
			% Read part of the modulation Segment
			fseek(fid, subSegmentOffset, 'bof');
			% Read first modulation period & duration of that period
			modTime = fread(fid, 1, 'double');
			modDur = fread(fid, 1, 'double');
			i = 1;
			while modDur(i, 1) ~ = -1
				i = i+1;
				fseek(fid, 34, 'cof');
				modTime(i, 1) = fread(fid, 1, 'double');
				modDur(i, 1) = fread(fid, 1, 'double');
			end
			if i > 1
				modDur(end) = round(time(end)/1000 - sum(modDur(1:end-1)), 0);
			else
				modDur = round(time(end)/1000, 0);
			end
		end

		% Find the CFilamentOffTimesStepd section of the header
		currSegmentOffset = 0;
		segmentFlag = 'CFilamentOffTimesStepd';
		flagData = bin2dec(dec2bin(segmentFlag));
		for i = 1:(offsets(3)-(offsets(2)+11))
			% Section starts with the string CFilamentOffTimesStepd
			if headerData(i:i+(length(segmentFlag)-1)) ==  flagData
				% Store the offset and stop search
				currSegmentOffset = offsets(2) + i + length(segmentFlag)-1;
				break;
			end
		end
		if (currSegmentOffset > 0)
			% Read part of the CFilamentOffTimesStepd section
			fseek(fid, currSegmentOffset + 1, 'bof');
			firstDelay = fread(fid, 1, 'double');
			% no seek, since last double drops us at the required offset
			filOff = fread(fid, 1, 'uint8');
			fseek(fid, 5, 'cof');
			nextDelay = fread(fid, 1, 'double');
		end

		% Correct time vector if acquisition was delayed
		if firstDelay > 0 && filOff ==  0 %if nextDelay is -1, then this is a simple delay, but if not then there is a table of filament on/off timings
			time(:) = time(:)+firstDelay*1000;
			acqDelay = firstDelay;
		else
			acqDelay = 0;
		end

	elseif magic ==  0
		errmsg = strcat(pegFile, ' corrupt - header data contains bad values.');
		error (errmsg)
	else
		errmsg = strcat(pegFile, ' does not appear to contain LECO .peg file data.');
		error (errmsg)
	end

	% Close the file
	fclose(fid);
end

% Move captured data into structure (memory heavy, but faster than direct struct referencing)
pegstruct.fileName = strrep(pegFile, "\\", "\");
pegstruct.acqDate = acqDate;
pegstruct.acqVer = acqVer;
pegstruct.dataRate = dataRate;
pegstruct.startMass = startMass;
for rng = 1:(massRange)
	pegstruct.masses(rng) = rng+(startMass)-1;
end
pegstruct.massRange = massRange;
pegstruct.acqDelay = acqDelay;
pegstruct.modTime = modTime;
pegstruct.modDur = modDur;
pegstruct.numScans = numScans;
pegstruct.time = time;
pegstruct.time1D = time/1000;
if procspec ==  1
	pegstruct.specdata = specdata;
	pegstruct.tic = tic;
	pegstruct.bpc = max(specdata,[],2);
end
pegstruct.msParams = msParams;
end
