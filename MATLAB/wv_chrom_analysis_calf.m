clear
clc
for z = 1:7
    for zz = 1:5
        testwvs=[688 795 808 830 860 915 976];
        if zz == 2
            addition = -10;
        elseif zz == 3
            addition = -5;
        elseif zz == 1
            addition = 0;
        elseif zz == 4
            addition = 5;
        elseif zz == 5;
            addition = 10;         
        end
        testwvs(z) = testwvs(z) + addition;
   
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
            
            clear fdpmfit
            % cd('C:\Users\hillb\Desktop\Better Data\PROCESSED')
            % [filename, pathname]=uigetfile('*MUA.asc','pick a measurement file')
            MUA_folder= dir('*_MUA_and_fit.asc');
            
            %my_loc=findstr(pathname,'\');
            %gimage.ID=pathname( (my_loc(end-1)+1):end-1);
            %gimage.rootbase = pathname(1: my_loc(end-1)-1 );
            %opt.datafile_dos={filename};
            
            % cd(pathname);
            filename = MUA_folder.name;
            data=dlmread(filename,'\t',6,0);
            fid=fopen(filename);
            C = textscan(fid, '%s','delimiter', '\n');
            testline = cell2mat(C{1,1}(1));
            if strcmp(testline(1:4),'wave')
                nameline = C{1}{1};
            else
            nameline=C{1}{6};
            end
            
            % Section for adding tags to data
            headers=textscan(nameline,'%s','delimiter','\t');
            clear DataTags
            for i = 1:1:length(headers{1,1})-3
                measureTag = cell2mat(headers{1,1}(i+1,1));
                if measureTag(1) == 'L'
                    DataTags(1,i) = 1;
                else
                    DataTags(1,i) = 0;
                end
                if measureTag(1) == 'R'
                    DataTags(2,i) = 1;
                else
                    DataTags(2,i) = 0;
                end
                DataTags(3,i) = str2num(measureTag(3:5));
                if measureTag(2) == 'N'
                    DataTags(3,i) = -DataTags(3,i);
                end
                DataTags(4,i) = str2num(measureTag(7:9));
                if measureTag(6) == 'N'
                    DataTags(4,i) = -DataTags(4,i);
                end
            end
            % Tag 1 is right Tag 2 is left Tag 3 is abdominal
            % Add quantification of spatial coordinate for exclusion or
            % inclusion of areas
            fclose(fid);
            
            wvs=data(:,1);
            
            [c,ia,indexes]=intersect(testwvs,wvs);
%             nfiles=(size(data,2))/2;
            nfiles=(size(data,2))-4; % -1 for wv's -2 for avg and std

            for j=1:nfiles
                muas=data(indexes,j+1);
                
                %% Physio Fit Settings
                physio.spec.baseline = 1;
                % bool
                % adds an additional constant to fit to during the physio fit
                % standard is on, we've found that this aids in fitting
                % chromophore values
                % note: this variable isn't used, was meant to replace spec.opt.baseline
                % to make code more clear
                
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
                % sets which chromophores from the physio.chrom.names array are used in the fittings
                
                physio.fittype = 3;
                % method for fitting physio
                %                 0   is simple least squares
                %                 1   is weighted least squares ( by fitted mua error)
                %                 2   simple SVD
                %                 3    constrained LSQ (positive only)
                %                 4    constrained LSQ
                
                
                
                fdpmchrom = physioSetup(physio, testwvs);
                [fdpmfit(j).phy, fdpmfit(j).fitmua, fdpmfit(j).ss] = physioFit(physio.fittype, fdpmchrom, muas', muas', 0);
            end
            
            %% bring in old fits
            clear broadband_fit_data_raw broadband_fit_data
            SUM_Folder = dir('*__SUM.asc');
            filename = SUM_Folder.name;
            fid=fopen(filename);
            D = textscan(fid, '%s', 'delimiter','\n');
            D_mat1 = cell2mat(D{1,1}(71,1));
            D_mat2 = cell2mat(D{1,1}(42,1));
            D_mat3 = cell2mat(D{1,1}(91,1));
            D_mat4 = cell2mat(D{1,1}(54,1));
            if (strfind(D_mat1, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw=dlmread(filename,'\t',[72 1 79 length(fdpmfit)]);
            elseif (strfind(D_mat2, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw =dlmread(filename,'\t',[43 1 50 length(fdpmfit)]);
            elseif (strfind(D_mat3, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw =dlmread(filename,'\t',[92 1 99 length(fdpmfit)]);
            elseif (strfind(D_mat4, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw =dlmread(filename,'\t',[56 1 63 length(fdpmfit)]);
            end
            
            % [filename, pathname]=uigetfile('*__SUM.asc','pick a measurement file')
           %if testdata=importdata(filename2)
%             dlmread(filename,'\t',[43 1 50 length(fdpm.fit)])
%             broadband_fit_data_raw=dlmread(filename,'\t',[72 1 79 length(fdpmfit)]);
            broadband_fit_data=broadband_fit_data_raw([1:4 6:8],:);
            fclose(fid);
            %% compile new fit data into matrix
            
            clear newFitArray oldFitArray
            newFitArray = zeros(length(fdpmchrom.names),nfiles);
            oldFitArray = zeros(length(fdpmchrom.names),nfiles);
            for j=1:nfiles
                for i =1:length(fdpmchrom.names)
                    newFitArray(i,j) = fdpmfit(j).phy(1,i);
                    oldFitArray(i,j) = broadband_fit_data(i,j);
                end               
            end
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