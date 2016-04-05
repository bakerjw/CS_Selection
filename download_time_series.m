function [  ] = download_time_series(outputDir, rec, Filename, dirLocation)
% Copy selected time series to the working directory, if they are available
% for direct download. For further details on file locations, see database
% metadata files.


% identify the relevant file format
if ~isempty(strfind(dirLocation{1}{1},'nga_files')) %% NGA-West1 data
    formatFlag = 1;
    % check for number of ground motion components
    if size(Filename,2) == 1 % only one ground motion component
        nComp = 1;
    else
        nComp = 2;
    end

elseif ~isempty(strfind(dirLocation{1}{1},'bbpvault')) %% BBP data
    formatFlag = 2;
    nComp = 1; % there is only one file for multiple components
    
elseif ~isempty(strfind(dirLocation{1}{1},'ngawest2')) %% NGA-West2 data
    display('The software currently cannot download NGA-West2 data');
    return
else
    display('Unknown database format');
    return
end

for i = 1:length(rec)
    for j = 1:nComp
        % Download the selected ground motion from the NGA database
        url = [char(dirLocation{rec(i)}) char(Filename{rec(i),j})]; % source data
        if nComp == 2
            filename = [outputDir '/GM' num2str(i) '_comp' num2str(j) '.txt']; % save file to the current directory
        else
            filename = [outputDir '/GM' num2str(i) '.txt']; % save file to the current directory
        end
        urlwrite(url,filename); % urlwrite(URL,filename) reads web content at the specified URL and saves it to the file specified by filename
    end
end

