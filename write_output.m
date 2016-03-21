function [  ] = write_output( rec, IMs, outputFile, getTimeSeries, Filename, dirLocation)
% Write a tab-delimited file with selected ground motions and scale factors. 
% For instructions on downloading the time histories, see the documentation
% files for each database. 

fin = fopen(outputFile,'w');

% print header information
fprintf(fin, '%s \n \n', getTimeSeries{1}, getTimeSeries{2}, getTimeSeries{3});
if size(Filename,2) == 1 % only one ground motion component
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \n','Record Number','Record Sequence Number','Scale Factor','File Name','URL');
else % two components
    fprintf(fin,'%s \t %s \t %s \t %s \t %s \t %s \t %s \n','Record Number','Record Sequence Number','Scale Factor','File Name Dir. 1','File Name Dir. 2', 'URL 1', 'URL 2');
end

% print record data
for i = 1 : length(rec)
    if size(Filename,2) == 1 
        fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \n',i,rec(i),IMs.scaleFac(i),char(Filename{rec(i)}),[char(dirLocation{rec(i)}) char(Filename{rec(i)})]); % Print relevant outputs
    else 
        fprintf(fin,'%d \t %d \t %6.2f \t %s \t %s \t %s \t %s \n',i,rec(i),IMs.scaleFac(i),char(Filename{rec(i),1}),char(Filename{rec(i),2}),[char(dirLocation{rec(i)}) char(Filename{rec(i),2})],[char(dirLocation{rec(i)}) char(Filename{rec(i),2})]);
    end
end

fclose(fin);


end

