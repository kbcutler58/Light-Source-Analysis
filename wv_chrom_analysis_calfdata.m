clear
clc
waves = 1;
for calf = 1:2
    
for waves = 1:29
    wvlengths = [...
        688	795	808	830	860	915	976;
        678	795	808	830	860	915	976;
        683	795	808	830	860	915	976;
        693	795	808	830	860	915	976;
        698	795	808	830	860	915	976;
        688	785	808	830	860	915	976;
        688	790	808	830	860	915	976;
        688	800	808	830	860	915	976;
        688	805	808	830	860	915	976;
        688	795	798	830	860	915	976;
        688	795	803	830	860	915	976;
        688	795	813	830	860	915	976;
        688	795	818	830	860	915	976;
        688	795	808	820	860	915	976;
        688	795	808	825	860	915	976;
        688	795	808	835	860	915	976;
        688	795	808	840	860	915	976;
        688	795	808	830	850	915	976;
        688	795	808	830	855	915	976;
        688	795	808	830	865	915	976;
        688	795	808	830	870	915	976;
        688	795	808	830	860	905	976;
        688	795	808	830	860	910	976;
        688	795	808	830	860	920	976;
        688	795	808	830	860	925	976;
        688	795	808	830	860	915	966;
        688	795	808	830	860	915	971;
        688	795	808	830	860	915	981;
        688	795	808	830	860	915	986];
    
    testwvs=wvlengths(waves,:);
    % testwvs = [688	795	808	830	860	915	976];
    cd 'C:\Users\Kyle\Downloads\Processed Data- 3 Diodes +SS'
    folders = dir('Processed Data- 3 Diodes +SS');
    cd('Processed Data- 3 Diodes +SS');
    %         folders = dir('PROCESSED_use these processed data');
    %         cd('PROCESSED_use these processed data')
    %         path1 = 'C:\Users\Kyle\Downloads\DOSI data for wv LG testing\PROCESSED_use these processed data';
    path1='C:\Users\Kyle\Downloads\Processed Data- 3 Diodes +SS\Processed Data- 3 Diodes +SS';
    for k = 1:length(folders)-2
        cd(strcat(path1,'\',folders(k+2).name));
        % cd((path(k,:)))%'C:\Users\Kyle\Downloads\DOSI data for wv LG testing\PROCESSED_use these processed data\6691-03_110803')
        clear fdpmfit fdpmfit2 physio PatientData
        
        
        MUA_folder= dir('*_MUA_and_fit.asc');
        if length(MUA_folder) > 1
            filename = MUA_folder(calf).name;
        else
            filename = MUA_folder.name;
        end
        
        
        % cd(pathname);
        %     filename = MUA_folder.name;
        data=dlmread(filename,'\t',6,0);
        fid=fopen(filename);
        C = textscan(fid, '%s','delimiter', '\n');
        testline = cell2mat(C{1,1}(1));
        if strcmp(testline(1:4),'wave')
            nameline = C{1}{1};
        else
            nameline=C{1}{6};
        end
        clear headers
        headers=textscan(nameline,'%s','delimiter','\t');
        accumulator = 0;
        clear header_index
        for i = 1:length(headers{1})-1
            measureTag = cell2mat(headers{1,1}(i+1,1));
            if ~strcmp(measureTag,'fit') && ~strcmp(measureTag,'average') && ~strcmp(measureTag,'standard deviation')
                header_index(i) = 1;
                if measureTag(1) == 'L'
                    DataTags(1,i-accumulator) = 1;
                else
                    DataTags(1,i-accumulator) = 0;
                end
                if measureTag(1) == 'R'
                    DataTags(2,i-accumulator) = 1;
                else
                    DataTags(2,i-accumulator) = 0;
                end
                DataTags(3,i-accumulator) = str2num(measureTag(5));
                %             if measureTag(2) == 'N'
                %                 DataTags(3,i-accumulator) = -DataTags(3,i-accumulator);
                %             end
                DataTags(4,i-accumulator) = str2num(measureTag(7));
                %             if measureTag(6) == 'N'
                %                 DataTags(4,i-accumulator) = -DataTags(4,i-accumulator);
                %             end
            else
                header_index(i) = 0;
                accumulator = accumulator + 1;
            end
        end
        
        % Tag 1 is left Tag 2 is right Tag 3 is xdist Tag 4 is ydist
        fclose(fid);
        wvs=data(:,1);
        [c,ia,indexes]=intersect(testwvs,wvs);
        nfiles=floor((size(data,2)-2)/2);
        data_fit = data(:,[false ~logical(header_index)]);
        data = data(:,[false logical(header_index)]);
        
        %% bring in old fits
        clear broadband_fit_data_raw broadband_fit_data
        SUM_Folder = dir('*__SUM.asc');
        if length(SUM_Folder) > 1
            filename = SUM_Folder(calf).name;
        else
            filename = SUM_Folder.name;
        end
        fid=fopen(filename);
        D = textscan(fid, '%s', 'delimiter','\n');
        D_mat3 = 0;
        D_mat1 = cell2mat(D{1,1}(71,1));
        D_mat2 = cell2mat(D{1,1}(42,1));
        D_mat4 = cell2mat(D{1,1}(54,1));
        D_mat5 = cell2mat(D{1,1}(36,1));
        if (length(D) > 91)
            D_mat3 = cell2mat(D{1,1}(91,1));
        end
        if (strfind(D_mat1, 'Conc (SSFDPM):')) == 1
            broadband_fit_data_raw=dlmread(filename,'\t',[72 1 79 nfiles]);
        elseif (strfind(D_mat2, 'Conc (SSFDPM):')) == 1
            broadband_fit_data_raw =dlmread(filename,'\t',[43 1 50 nfiles]);
        elseif (strfind(D_mat3, 'Conc (SSFDPM):')) == 1
            broadband_fit_data_raw =dlmread(filename,'\t',[92 1 99 nfiles]);
        elseif (strfind(D_mat4, 'Conc (SSFDPM):')) == 1
            broadband_fit_data_raw =dlmread(filename,'\t',[56 1 63 nfiles]);
        elseif (strfind(D_mat5, 'Conc (SSFDPM):')) == 1
            broadband_fit_data_raw =dlmread(filename,'\t',[37 1 44 nfiles]);
        end
        
        % [filename, pathname]=uigetfile('*__SUM.asc','pick a measurement file')
        %if testdata=importdata(filename2)
        %             dlmread(filename,'\t',[43 1 50 length(fdpm.fit)])
        %             broadband_fit_data_raw=dlmread(filename,'\t',[72 1 79 length(fdpmfit)]);
        %             broadband_fit_data=broadband_fit_data_raw([1:4 6:8],:);
        broadband_fit_data=broadband_fit_data_raw([1:4 6:8],:);
        fclose(fid);
        %%
        for j=1:nfiles
            %% Main Function
            
            muas=data(indexes,j);
            physio.spec.baseline = 1;
            physio.chrom.file = 'chromophores_Zijlstra.600.txt';
            % files located in ssfdpmPro\chromophore_files
            % contain extiction coefficients of chromophores
            % different files can be used for different chromophores
            % can also have extiction coefficients of same chromophores found from
            % different sources
            physio.chrom.names = {'HbO2', 'Hb', 'h2oFrac', 'fatFrac', 'Met','CO','Evans','MB','bkgd'};
            % names of chromophores found in chromophore file
            physio.chrom.units = {'uM', 'uM', '%', '%', 'uM','uM','uM','uM','uM'};
            % units to use for chomophores found in chomophore file
            physio.chrom.mult = [10^3; 10^3; 100; 100; 10^3;10^3;10^3;10^3;10^3];
            % conversion multiplier to use after fitting chromophores concentrations
            physio.chrom.mins = [0;0;0;0;0;0;0;0];
            physio.chrom.maxs = [.1;.1;1;1;.1;.1;0;0];
            % if physio.fittype is constrained (see below) these variables set the
            % lower and upper bounds of the fit
            physio.chrom.selected = [1,1,1,1,0,0,0,0];
            physio.fittype = 3;
            
            fdpmchrom = physioSetup(physio, testwvs,0);
            [fdpmfit(j).phy, fdpmfit(j).fitmua, fdpmfit(j).ss] = physioFit(physio.fittype, fdpmchrom, muas', muas', 0);
            %                 physio.chrom.maxs = [.1;.1;1;1;.1;.1;0;0];
            %                 muas=data(indexes,j*2);
            fdpmchrom2 = physioSetup(physio, wvs, 0);
            concentrations = fdpmfit(j).phy(1,1:4);
            [~, fdpmfit2(j).fitmua] = physioFit3(physio.fittype, fdpmchrom2, muas', muas', 0, concentrations);
            mua_fit_spec = fdpmfit2(j).fitmua;
            
            PointLabel = strcat('Point',num2str(j,'%03d'));
            PatientData.(PointLabel).wavelengths = testwvs;
            PatientData.(PointLabel).muaSpectrum(1,:)=wvs';
            PatientData.(PointLabel).muaSpectrum(2,:)=mua_fit_spec;
            PatientData.(PointLabel).muaSpectrum(3,:)=data(:,j);
            PatientData.(PointLabel).muaSpectrum(4,:)=data_fit(:,j);
            PatientData.(PointLabel).pointlabel = headers{1,1}([false logical(header_index)]);
            PatientData.(PointLabel).tag = DataTags(:,j);
            PatientData.(PointLabel).chromo(1,:) = fdpmfit(1,j).phy(1,:);
            PatientData.(PointLabel).chromo(2,:) = broadband_fit_data(:,j)';
            PatientData.(PointLabel).chromo(3,:) = broadband_fit_data(:,j)'-fdpmfit(1,j).phy(1,:);
            PatientData.(PointLabel).chromo(4,:) = ((broadband_fit_data(:,j)'-fdpmfit(1,j).phy(1,:))./broadband_fit_data(:,j)')*100;
        end
        if calf == 1
            OutputFile = strcat(folders(k+2).name,'_L','_waves',num2str(waves,'%02d'),'.mat');
        else
            OutputFile = strcat(folders(k+2).name,'_R','_waves',num2str(waves,'%02d'),'.mat');
        end
        cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Calf Results');
        save(OutputFile,'PatientData')
    end
    
end
end
%{
%% compile new fit data into matrix
    
    %% db files
    % outputfile = strcat('wvanalysis_SUM_itr',num2str(((z-1)*5)+zz),'.asc');
    cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
    fileoutputname = strcat('Patient',filename(1:16),'_ABS3.asc');
    fid=fopen(fileoutputname,'a');
    % fid=fopen(outputfile,'w');
    % fid=fopen(outputfile,'a');
    %fprintf(fid, 'patientID \tdate \tposition \twv \tmua \tdmua \tmus \tdmus \tfitmethod \tcomment\n');
    %for i = 1:length(chrom.names),disp(sprintf('[%s]\t= %3.2f +/- %2.3f %s', chrom.names{i}, conc(i), conc_err(i), chrom.units{i})); end
    fprintf(fid,'\t%s',num2str(testwvs));
    fprintf(fid,'\n');
    for j=1:nfiles
        fprintf(fid,'\t%s',headers{1}{j+1});
    end
    fprintf(fid,'\n');
    for t = 1:4
        for j=1:nfiles
            fprintf(fid,'\t%f',DataTags(t,j));
        end
        fprintf(fid,'\n');
    end
    for i=1:length(fdpmchrom.names)
        fprintf(fid,'%s',fdpmchrom.names{i});
        for j=1:nfiles
            fprintf(fid,'\t%f',fdpmfit(j).phy(1,i));
        end
        fprintf(fid,'\n');
        for j=1:nfiles
            fprintf(fid,'\t%f',broadband_fit_data(i,j));%fdpmfit(j).phy(1,i));
        end
        fprintf(fid,'\n');
        for j=1:nfiles
            fprintf(fid,'\t%f',(fdpmfit(j).phy(1,i))-(broadband_fit_data(i,j)));%fdpmfit(j).phy(1,i));
        end
        fprintf(fid,'\n');
        for j=1:nfiles
            fprintf(fid,'\t%f',(((fdpmfit(j).phy(1,i))-(broadband_fit_data(i,j)))/(broadband_fit_data(i,j))*100));%fdpmfit(j).phy(1,i));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    
    fclose(fid);
    outputfile = strcat('wvanalysis_Feb24_1_results_ABS3','.asc');
    fid = fopen(outputfile,'a');
    fprintf(fid,'\t%s', folders(k+2).name);
    fprintf(fid,'\t%s', num2str(testwvs));
    %             fprintf(fid,'\t%s', 'Analysis Total');
    for i=1:length(fdpmchrom.names)
        fprintf(fid,'\t%f', mean(newFitArray(i,:)-oldFitArray(i,:)));
        fprintf(fid,'\t%f', std(newFitArray(i,:)-oldFitArray(i,:)));
        fprintf(fid,'\t%f', mean(((newFitArray(i,:)-oldFitArray(i,:))./oldFitArray(i,:))*100));
        fprintf(fid,'\t%f', std(((newFitArray(i,:)-oldFitArray(i,:))./oldFitArray(i,:))*100));
    end
    fprintf(fid,'\n');
    fclose(fid);
end

end
end
%}