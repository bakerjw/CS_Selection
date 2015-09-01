% Please run this script after "Select_Ground_Motions" to format outputs
% (modified for CEE 385 interface with IIIDAP v1.5)

% Ting Lin 
% 10/9/2012

for i = 1 : length(finalRecords)
    % Download the selected ground motion from the NGA database
    rec = finalRecords(i);
    url = ['http://peer.berkeley.edu/nga_files/ath/' Filename{rec}(1:end-3) 'AT2']; % URL of selected ground motion
    filename = ['GM' num2str(i) '.AT2']; % save file to the current directory
    urlwrite(url,filename); % urlwrite(URL,filename) reads web content at the specified URL and saves it to the file specified by filename
    
    % Read in relevant data
    file_in = fopen(filename, 'r');    % open the file for reading
    for li=1:3
        ans = fgetl(file_in);   % discard the first 3 lines
    end
    [ans ans] = fscanf(file_in, '%5c', 1);   
    [record_dt ans] = fscanf(file_in, '%f', 1);  % read in time step, dt
    ans = fgetl(file_in);               % read in the rest of line 4
    [A, np] = fscanf(file_in, '%f');    % A is a vector of all accelerations, np is length of A = number of points
    fclose(file_in);                    % close files
    
    % Store time increments and number of points
    GM_dt(i,1) = record_dt;
    GM_npts(i,1) = np;
    
    % Save acceleration time history for each selected ground motion
    fname=['GroundMotion' num2str(i) '.th'];
    save (fname, 'A', '-ASCII')
end

% Save time steps for all selected ground motions
fname=['GroundMotionTimeStep.out'];
save (fname, 'GM_dt', '-ASCII')
% Save number of points for all selected ground motions
fname=['GroundMotionTotalPoints.out'];
save (fname, 'GM_npts', '-ASCII')