function [  ] = write_output(recIdx, recNum, IMs, outputDir, outputFile, dirLocation, metadata)
% Write a tab-delimited file with selected ground motions and scale factors. 
% For instructions on downloading the time histories, see the documentation
% files for each database. 

% Create directory for outputs, if it doesn't yet exist
if ~exist(outputDir, 'dir') 
     mkdir(outputDir)
end

fin = fopen([outputDir '/' outputFile],'w');

% print header information
fprintf(fin, '%s \n \n', metadata.getTimeSeries{1}, metadata.getTimeSeries{2}, metadata.getTimeSeries{3});
if size(metadata.Filename,2) == 1 % only one ground motion component
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \t %s \n','Record Number','Record Sequence Number','Scale Factor','Component Number','File Name','URL');
else % two components
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \t %s \t %s \n','Record Number','Record Sequence Number','Scale Factor','File Name Dir. 1','File Name Dir. 2', 'URL 1', 'URL 2');
end

% print record data
for i = 1 : length(recIdx)
    if size(metadata.Filename,2) == 1 % only one ground motion component
        fprintf(fin,'%d \t %d \t %6.2f \t %d \t %s \t %s \n',i,recNum(recIdx(i)),IMs.scaleFac(i), metadata.compNum(recIdx(i)), char(metadata.Filename{recIdx(i)}),[char(dirLocation{recIdx(i)}) char(metadata.Filename{recIdx(i)})]); % Print relevant outputs
    else 
        fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \t %s \t %s \n',i,recIdx(i),IMs.scaleFac(i),char(metadata.Filename{recIdx(i),1}),char(metadata.Filename{recIdx(i),2}),[char(dirLocation{recIdx(i)}) char(metadata.Filename{recIdx(i),2})],[char(dirLocation{recIdx(i)}) char(metadata.Filename{recIdx(i),2})]);
    end
end

fclose(fin);


end

