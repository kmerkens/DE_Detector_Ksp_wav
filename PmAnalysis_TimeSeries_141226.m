%%PmAnalysis_TimeSeries_XXXXXX.m
%%Code to calculate the running mean of detections for 14 day periods for
%%each deployment at each site.  
%%Uses SampleData output from PmAnalysis_140701.m, which contains the 
%%matlab datenum for every 5 minute cycle with detections per site. 
%%Also uses the effort data from PIR_HARP_data_summary.xls to identify the
%%start and end of each deployment to correctly calculate the running
%%means. 
%%To change the length of the smooth, adjust parameter ss, line 174


clearvars  %remove variables from whatever was run before
close all

%Make cell array for storing all data together for future plotting
allsitesdata = {};


%% Identify the output files generated for all of the regions, all saved in
%#SampleData directory
dirname = ('C:\Users\Karlina.Merkens\Documents\SpermWhales\MATLAB_output\#SampleData');
dirdata = what(dirname);
matfiles = dirdata.mat;
nummatfiles = size(matfiles,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get data sheet
%user input version
% % [effsheetname, effsheetpath] = uigetfile('C:\Users\Karlina.Merkens\Documents\PIR_Misc\*.xls',...
% %     'Select the file that contains the effort information');
% % efffile = [effsheetpath, effsheetname];
%hardcoded version
efffile = 'C:\Users\Karlina.Merkens\Documents\PIR_Misc\PIR_HARP_data_summary.xls';

[effnums, effstrings, effraw] = xlsread(efffile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract deployment names (taken from dutycycleanalysis_140506.m, see there
%for annotation)
reg_dep_sites = char(effstrings(:,2));  %Character array of SIO names, 
regions = upper(reg_dep_sites(:,1:3));
numdeps = length(reg_dep_sites); %how many to loop through
deplymts = NaN(numdeps,1); %initialize
for dp = 1:1:numdeps 
    trim_reg_dep_site = strtrim(reg_dep_sites(dp,:)); %trim off white space
    if isempty(trim_reg_dep_site)
        %display('Skipping row - no deployment data')
        continue
    end
    namesplit = textscan(trim_reg_dep_site,'%s','delimiter','-'); %divides 
    %the name at the - that falls between the deployment number and the
    %site name.
    regiondep = namesplit{1,1}{1,1};
    depnumber = str2double(regiondep(end-1:end));
    deplymts(dp,1) = depnumber;
    sites(dp,1) = trim_reg_dep_site(end);
end
sites = upper(sites);


%% For each matfile, proceeed with analysis
finalindex = 1;

for mf = 1:nummatfiles
    matfile = char(matfiles(mf));
    mattype = matfile(18:20); %to identify whether twosample ddata or subsampled data.
    matfileopen = [dirname,'\',matfile];
    load(matfileopen)
    
    subsamp = strfind(matfile,'Sub');
    if ~isempty(subsamp)
        continue %move on to another file - this is not two sample.
    end
    
    region = upper(matfile(12:14)); %3 characters of region
    site = matfile(16); %1 character, last letter of the site code
    site = upper(site);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Find the matching rows of effort data
    regionmatch = ismember(regions, region(1:3), 'rows');
    sitematch = ismember(sites, site, 'rows');
    matches = regionmatch + sitematch; %add across the matched vectors 
    %and those that equal 2 are the matches
    match = find(matches == 2);
    howmany = size(match,1); %just check that there is at least one that matches
    %Check that the duty cycle isn't too long (sai-A-01)
   prune = [];
    for h = 1:howmany
         dutycycle = effraw{match(h),14};
        if dutycycle > 30
            prune = [prune, h];
            display(['Removing one because duty cycle is too long',char(effraw(match(h)))])
        end
    end
    match(prune) = [];
    %check again
    howmany = size(match,1); 
    if howmany < 1
        display('Uh Oh! There was no match for this region & site in the database')
        return
    end

    deprow = match; %these are the rows that are the deployments of interest

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract start/end times of each deployment
    %Start dates are in column 21, ends in 22. Use "raw" because the "nums" 
    %chops of rows without any numbers, which could get confusing, whereas 
    %raw stores everything.
    %Then convert to matlab date number
    numregdeps = size(deprow,1);
    startdate = zeros(numregdeps,1);
    startdatemat = zeros(numregdeps,1);
    enddate = zeros(numregdeps,1);
    enddatemat = zeros(numregdeps,1);
    i = 1;
    for dr = 1:1:numregdeps
        if ~isnan(effraw{match(dr), 21})
            startdate(i,1) = effraw{match(dr), 21};
            startdatemat(i,1) = startdate(i,1) + 693960;
            enddate(i,1) = effraw{match(dr), 22};
            enddatemat(i,1) = enddate(i,1) + 693960; 
            i = i + 1;
        else
            display(['No start date this deployment: ',  effraw{match(dr), 2}]);
            startdate(i) = [];
            startdatemat(i) = [];
            enddate(i) = [];
            enddatemat(i) = [];
            continue
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Fill in list of dates, to identify the gaps between deployments
    %alldates = [];
    startdatesround = ceil(startdatemat);
    enddatesround = floor(enddatemat);
    alldates = [min(startdatesround):1:max(enddatesround)]';
    
    effdates = [];
    numeffdeps = size(startdatemat,1);
    for nd = 1:1:numeffdeps
        depstart = ceil(startdatemat(nd)); %round up in case not full day
        depend = floor(enddatemat(nd)); %round down in case not full day
        depalldates = [depstart:1:depend]';
        effdates = [effdates; depalldates];
    end

    %Sort effalldates, just in case they're not chronological
    effdates = sort(effdates);
    %Find any dates in the long time series that aren't in the deployments
    %First compare the two vectors
    gaplogic = ismember(alldates,effdates);
    %Find where there are 0s
    gaprows = find(0 == gaplogic(:,1));
    %Make vector of 0s for second row of data matrix
    numalldates = size(alldates,1);
    alldateszeros = zeros(numalldates,1);
    dataperdate = [alldates, alldateszeros];
    dataperdate(gaprows,2) = NaN;
    
    dailyeffort = (60*24)/30; %number of 5 minute cycles per day, for two 
    %sample - 30 minute - data. 
    
    numcells = size(dataarray,2);
    for sc = 1:1:numcells
        ncycles = size(dataarray{1,sc},1);
        for nc = 1:1:ncycles %For each cycle recorded with detections...
            cycledate = floor(dataarray{1,sc}(nc,1)); %round down to the day
            daterow = find(cycledate == dataperdate); %Find the row in dataperdate
            %that matches with that rounded date
            dataperdate(daterow,2) = dataperdate(daterow,2)+1; %add one to that date
        end
    end

    dataperdate(:,3) = dataperdate(:,2)/dailyeffort; %divde by the total cycles 
    %possible to get the proportion.
    
    %Add padding NaNs before and after the data, to make running nanmean
    %work well. 
    %ss is half the length of the smooth minus 1 (so equal amounts on each side)
    ss = 30;
    smoothstr = num2str(ss*2+1);
    display(['Smooth length is: ',num2str(smoothstr)])
    padnans = nan(ss,1);
    beforepadstart = min(alldates) - ss;
    beforepadend = min(alldates)-1;
    beforedates = [beforepadstart:1:beforepadend]';
    beforepad = [beforedates,padnans,padnans];
    afterpadstart = max(alldates) + 1;
    afterpadend = max(alldates) + ss;
    afterpaddates = [afterpadstart:1:afterpadend]';
    afterpad = [afterpaddates,padnans,padnans];
    dataperdatepadded = [beforepad; dataperdate; afterpad];
    
    %Now loop through each date, and calculate a nanmean for the
    %surrounding ss*2 days. 
    numdates = size(dataperdate,1);
    ssfull = (ss*2)+1;
    for pd = 1:1:(numdates) 
        endpoint = pd + ssfull - 1;
        chunk = dataperdatepadded(pd:endpoint,3);
        datenanmean = nanmean(chunk);
        dataperdate(pd,4) = datenanmean;
    end
    
    %For triangular smooth (two passes) repeat
    dataperdatepadded(ss+1:end-ss,4) = dataperdate(:,4);
    for pd = 1:1:(numdates) 
        endpoint = pd + ssfull - 1;
        chunk = dataperdatepadded(pd:endpoint,4);
        datenanmean = nanmean(chunk);
        dataperdate(pd,5) = datenanmean;
    end
    
% % %     %Or try the fastsmooth.m function
% % %     SmoothY = fastsmooth(dataperdate(:,3),ssfull,1,1);
% % %     dataperdate(:,5) = SmoothY;
% % %     
    %paste nans over the gaps again, to remove any nanmeans that might have
    %spread into the gap
    dataperdate(gaprows,4:5) = NaN;
    
    %%%Save the data for later plotting.
    allsitesdata{finalindex,1}.name = [region,'-',site];
    allsitesdata{finalindex,1}.data = dataperdate;
    finalindex = finalindex + 1;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make some plots, just for fun
%These two have them all in chronological order
for fignum = 1:3
    legendnames = {};
    %cc = jet(13);
    %For Color
    whitebg('w'); %make sure background color is set standard
    cc = [0,1,1; 0,0,0; 1,0,1; 0,0,0; 0,1,1; .75,.5,0; 0.5,1,.5; .5,1,.5; ...
        1,1,0; 0,0,1; 1,.75,0; .3,.6,.3; 1,0,0];
    lw = [3;3;3;3;3;3;3;3;3;3;3;3;3];
    ll = ['-'; '-' ;'-'; '-'; '-'; '-'; '-'; '-'; '-'; '-'; '-'; '-'; '-'];
    %For B&W
%     cc = [.8,.8,.8; 0,0,0; 0,0,0; 0,0,0; .5,.5,.5; .5,.5,.5; .2,.2,.2; .8,.8,.8; ...
%         .5,.5,.5; .8,.8,.8; 0,0,0; .8,.8,.8; 0,0,0];
%     lw = [3;3;3;2;2;2;2;2;2;3;2;3;3];
%     ll = ['--'; '--' ;'- '; '- '; '- '; '--'; '--'; '- '; '- '; '--'; '- '; '- '; '--'];
%     
    
    figure(1)
    for jj = 1:finalindex-1
       xdata = allsitesdata{jj,1}.data(:,1);
       ydata = allsitesdata{jj,1}.data(:,4);
       pp(jj,:) = plot(xdata, ydata,'color',cc(jj,:), 'LineStyle', ll(jj,:),...
           'LineWidth', lw(jj,:));
       hold on
       dataname = allsitesdata{jj,1}.name;
       legendnames{jj} = [dataname];
    end
    %title('''Boxcar'' sliding-average smooth (one pass)');
    ylabel({'Mean proportion of cycles per day with detections:';[smoothstr,' day smooth']});
    
    
    if fignum == 1
        %ylim([-0.005 0.06]);
        ylim([-0.005 0.301]);
        current_day = datenum('Jan 01 2005'); %the date the plot starts
        [Y, M, D, H, MN, S] = datevec(current_day);
        current_day = addtodate(current_day, -D + 1,'day');
        last_day = datenum('Dec 31 2007');
        legend([pp(1,:), pp(3,:), pp(7,:), pp(9,:)],...
            'CSM','HAW', 'LSM-S', 'PAL-WT',...
            'Location','northwest')
    elseif fignum == 2
        ylim([-0.005 0.301]);
        current_day = datenum('Jan 1 2008'); %the date the plot starts
        [Y, M, D, H, MN, S] = datevec(current_day);
        current_day = addtodate(current_day, -D + 1,'day');
        last_day = datenum('Dec 31 2010');
        legend([pp(3,:), pp(4,:), pp(6,:), pp(8,:), pp(9,:), pp(10,:), pp(13,:)],...
            'HAW','KAU','LSM-D','PAL-NS','PAL-WT','PHR','WAK',...
            'Location','northwest')
    elseif fignum == 3
        ylim([-0.005 0.301]);
        current_day = datenum('Jan 1 2011'); %the date the plot starts
        [Y, M, D, H, MN, S] = datevec(current_day);
        current_day = addtodate(current_day, -D + 1,'day');
        last_day = datenum('Dec 31 2013');
        legend([pp(2,:), pp(3,:),  pp(5,:), pp(10,:), pp(11,:), pp(12,:), pp(13,:)],...
            'EQU','HAW','KIN', 'PHR', 'SAI', 'TIN','WAK')
    end
        
    
    xlim([current_day last_day]);
    xtick = [current_day];
    xstep_length = 12; %one tick per year
    xstep_unit = 'month';
    while (last_day > current_day)
        xtick = [xtick addtodate(current_day, xstep_length, xstep_unit)];
        current_day = addtodate(current_day, xstep_length, xstep_unit);
    end

    ts = figure(1);
    ts = gca;
    set(ts,'XTick', xtick,'XTickLabel', datestr(xtick,'yyyy')); 
    %,'FontSize',14
    set(gca, 'Ticklength', [0 0])
    plot(get(gca,'xlim'), [0 0], ':'); % adds horizontal line at 0. Adapts to x limits of current axes
    set(gcf, 'renderer', 'zbuffer'); %removes the exponent from the date
    %set size
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 10 3.5]); %x_width=10in y_width=3.5in
    whitebg %changes background to black to make colors more visible

    datestamp = datestr(now,'yymmdd');
    if fignum == 1
        yearsname = '2005-2007';
    elseif fignum == 2
        yearsname = '2007-2010';
    elseif fignum == 3
        yearsname = '2011-2013';
    end
    savefigt = (['TimeSeries_chronological_boxcar',yearsname,'-',datestamp,'.tiff']);
    %savefigj = (['TimeSeries_chronological_boxcar',datestamp,'.jpg']);
    plotsavedir = ('C:\Users\Karlina.Merkens\Documents\SpermWhales\MATLAB_output\#TimeSeries_Plots');
    savepath = ([plotsavedir,'\',savefigt]);
    set(gcf, 'InvertHardCopy', 'off');
    saveas(ts,savepath);
%     savepath = ([plotsavedir,'\',savefigj]);
%     saveas(ts,savepath);
    close all
end

% %%%
% %Make a plot of the smooth data
% whitebg
% figure(2)
% for kk = 1:finalindex-1
%    xdata = allsitesdata{kk,1}.data(:,1);
%    ydata = allsitesdata{kk,1}.data(:,5);
%    plot(xdata, ydata,'color',cc(kk,:));
%    hold on
% end
% title('''Triangular'' smooth (two passes)');
% ylabel(['Mean proportion of cycles per day with detections: ',smoothstr,' day smooth']);
% ylim([-0.05 0.3]);
% legend(legendnames)
% 
% current_day = datenum('Jan 01 2005'); %the date the plot starts
% [Y, M, D, H, MN, S] = datevec(current_day);
% current_day = addtodate(current_day, -D + 1,'day');
% last_day = datenum('Nov 1 2013');
% xlim([current_day last_day]);
% xtick = [current_day];
% xstep_length = 12; %one tick per year
% xstep_unit = 'month';
% while (last_day > current_day)
%     xtick = [xtick addtodate(current_day, xstep_length, xstep_unit)];
%     current_day = addtodate(current_day, xstep_length, xstep_unit);
% end
% 
% ts = figure(2);
% ts = gca;
% set(ts,'XTick', xtick,'XTickLabel', datestr(xtick,'yyyy')); 
% %,'FontSize',14
% set(gca, 'Ticklength', [0 0])
% plot(get(gca,'xlim'), [0 0], ':'); % adds horizontal line at 0. Adapts to x limits of current axes
% set(gcf, 'renderer', 'zbuffer'); %removes the exponent from the date
% %set size
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 10 3.5]); %x_width=10in y_width=3.5in
% whitebg
% 
% savefigt = (['TimeSeries_chronological_triangle',datestamp,'.tiff']);
% savefigj = (['TimeSeries_chronological_triangle',datestamp,'.jpg']);
% plotsavedir = ('C:\Users\Karlina.Merkens\Documents\SpermWhales\MATLAB_output\#TimeSeries_Plots');
% savepath = ([plotsavedir,'\',savefigt]);
% set(gcf, 'InvertHardCopy', 'off');
% saveas(ts,savepath);
% savepath = ([plotsavedir,'\',savefigj]);
% saveas(ts,savepath);
% close all
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %These two have them all all starting at 0
% whitebg
% figure(3)
% for jj = 1:finalindex-1
%    ydata = allsitesdata{jj,1}.data(:,4);
%    plot(ydata,'color',cc(jj,:));
%    hold on
% end
% title('''Boxcar'' sliding-average smooth (one pass)');
% xlabel('Days from start of first deployment');
% ylabel(['Mean proportion of cycles per day with detections: ',smoothstr,' day smooth']);
% ylim([-0.05 0.3]);
% legend(legendnames)
% 
% ts = figure(3);
% ts = gca;
% set(gca, 'Ticklength', [0 0])
% set(gcf, 'renderer', 'zbuffer'); %removes the exponent from the date
% %set size
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 10 3.5]); %x_width=10in y_width=3.5in
% whitebg
% 
% savefigt = (['TimeSeries_from0_boxcar',datestamp,'.tiff']);
% savefigj = (['TimeSeries_from0_boxcar',datestamp,'.jpg']);
% plotsavedir = ('C:\Users\Karlina.Merkens\Documents\SpermWhales\MATLAB_output\#TimeSeries_Plots');
% savepath = ([plotsavedir,'\',savefigt]);
% set(gcf, 'InvertHardCopy', 'off');
% saveas(ts,savepath);
% savepath = ([plotsavedir,'\',savefigj]);
% saveas(ts,savepath);
% close all
% 
% %Make a plot of the smooth data
% whitebg
% figure(4)
% for kk = 1:finalindex-1
%    ydata = allsitesdata{kk,1}.data(:,5);
%    plot(ydata,'color',cc(kk,:));
%    hold on
% end
% %title('Smoothed using fastsmooth.m')
% title('''Triangular'' smooth (two passes)');
% xlabel('Days from start of first deployment');
% ylabel(['Mean proportion of cycles per day with detections: ',smoothstr,' day smooth']);
% ylim([-0.05 0.3]);
% legend(legendnames)
% 
% ts = figure(4);
% ts = gca;
% set(gca, 'Ticklength', [0 0])
% set(gcf, 'renderer', 'zbuffer'); %removes the exponent from the date
% %set size
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 10 3.5]); %x_width=10in y_width=3.5in
% whitebg
% 
% savefigt = (['TimeSeries_from0_triangle',datestamp,'.tiff']);
% savefigj = (['TimeSeries_from0_triangle',datestamp,'.jpg']);
% plotsavedir = ('C:\Users\Karlina.Merkens\Documents\SpermWhales\MATLAB_output\#TimeSeries_Plots');
% savepath = ([plotsavedir,'\',savefigt]);
% set(gcf, 'InvertHardCopy', 'off');
% saveas(ts,savepath);
% savepath = ([plotsavedir,'\',savefigj]);
% saveas(ts,savepath);
% close all
% 
% 
% 
% 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %Save allsitesdata as .mat
% % % % savefilename = (['DataPerDate_',smoothstr,'-daysmooth_',datestamp]);
% % % % savedirname = ('C:\Users\Karlina.Merkens\Documents\SpermWhales\MATLAB_output\#DataPerDate');
% % % % savefilename_mat = [savedirname,'\', savefilename,'.mat'];
% % % % if exist(savefilename_mat) == 2
% % % %     delete(savefilename_mat)
% % % % end
% % % % save(savefilename_mat, 'allsitesdata','-mat');
% % % % 
% % % % %To save as xls, Convert to xls datenum, save each as a new sheet
% % % % savefilename = (['DataPerDate_',smoothstr,'-daysmooth_',datestamp]);
% % % % savefilename_xls = [savedirname,'\', savefilename,'.xls'];
% % % % if exist(savefilename_xls) == 2
% % % %     delete(savefilename_xls)
% % % % end
% % % % 
% % % % for sh = 1:finalindex-1
% % % %     allsitesdata_xls =  allsitesdata{sh,1}.data(:,1) - 693960;
% % % %     allsitesdata_xls(:,2:4) = allsitesdata{sh,1}.data(:,2:4);
% % % %     dataname = allsitesdata{sh,1}.name;
% % % %     xlswrite(savefilename_xls, allsitesdata_xls, dataname)
% % % % end
% 

