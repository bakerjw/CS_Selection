function [  ] = download_time_series(outputDir, rec, metadata)
% Copy selected time series to the working directory, if they are available
% for direct download. For further details on file locations, see database
% metadata files.


% identify the relevant file format
if ~isempty(strfind(metadata.dirLocation{1},'nga_files')) %% NGA-West1 data
    % check for number of ground motion components
    if size(metadata.Filename,2) == 1 % only one ground motion component
        nComp = 1;
    else
        nComp = 2;
    end
    
    for i = 1:length(rec)
        for j = 1:nComp
            % Download the selected ground motion from the NGA database
            url = [char(metadata.dirLocation{rec(i)}) char(metadata.Filename{rec(i),j})]; % source data
            if nComp == 2
                filename = [outputDir '/GM' num2str(i) '_comp' num2str(j) '.txt']; % save file to the current directory
            else
                filename = [outputDir '/GM' num2str(i) '.txt']; % save file to the current directory
            end
            urlwrite(url,filename); % urlwrite(URL,filename) reads web content at the specified URL and saves it to the file specified by filename
        end
    end

elseif ~isempty(strfind(metadata.dirLocation{1},'bbpvault')) %% BBP data
    display('BroadBand Platform data are currently not available for download');

elseif ~isempty(strfind(metadata.dirLocation{1},'ngawest2')) %% NGA-West2 data
    display('NGA-West2 data are currently not available for download');
    return
elseif ~isempty(strfind(metadata.dirLocation{1},'CyberShake')) %% CyberShake data
    if size(metadata.Filename,2) == 1 % only one ground motion component
        nComp = 1;
    else
        nComp = 2;
    end
    
    for i = 1:length(rec)
        for j = 1:nComp
            % Download the selected ground motion from the NGA database
            sourceFile = [char(metadata.dirLocation{rec(i)}) char(metadata.Filename{rec(i),j})]; % source data
            if nComp == 2
                filename = [outputDir '/GM' num2str(i) '_comp' num2str(j) '.txt']; % save file to the current directory
            else
                filename = [outputDir '/GM' num2str(i) '.txt']; % save file to the current directory
            end
            copyfile(sourceFile,filename); % urlwrite(URL,filename) reads web content at the specified URL and saves it to the file specified by filename
        end
    end
   
    
else
    display('Unknown database format');
    return
end



