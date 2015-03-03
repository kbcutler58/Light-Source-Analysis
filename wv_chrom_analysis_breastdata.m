clear
clc
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

   
        cd 'C:\Users\Kyle\Downloads\DOSI data for wv LG testing'
        folders = dir('PROCESSED_use these processed data');
        cd('PROCESSED_use these processed data');
%         folders = dir('PROCESSED_use these processed data');
%         cd('PROCESSED_use these processed data')
%         path1 = 'C:\Users\Kyle\Downloads\DOSI data for wv LG testing\PROCESSED_use these processed data';
        path1='C:\Users\Kyle\Downloads\DOSI data for wv LG testing\PROCESSED_use these processed data';
        for k = 1:length(folders)-2
            
            %% Load in Mua Spectrum and Previous Concentration Values
            cd(strcat(path1,'\',folders(k+2).name));
            % cd((path(k,:)))%'C:\Users\Kyle\Downloads\DOSI data for wv LG testing\PROCESSED_use these processed data\6691-03_110803')
            clear fdpmfit fdpmfit2 physio PatientData
            % Find Folder with _MUA.asc
            MUA_folder= dir('*_MUA.asc');
            filename = MUA_folder.name;
            clear data
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
            data_fitting = 0;
            for i = 1:length(headers{1})-1
                measureTag = cell2mat(headers{1,1}(i+1,1));
                if strcmp(measureTag,'fit')
                    data_fitting = 1;
                end
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
                    DataTags(3,i-accumulator) = str2num(measureTag(3:5));
                    if measureTag(2) == 'N'
                        DataTags(3,i-accumulator) = -DataTags(3,i-accumulator);
                    end
                    DataTags(4,i-accumulator) = str2num(measureTag(7:9));
                    if measureTag(6) == 'N'
                        DataTags(4,i-accumulator) = -DataTags(4,i-accumulator);
                    end
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
            if data_fitting == 1;
            data_fit = data(:,[false ~logical(header_index)]);
            end
            data = data(:,[false logical(header_index)]);
%             PatientData = struct();
            
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
                broadband_fit_data_raw=dlmread(filename,'\t',[72 1 79 nfiles]);
            elseif (strfind(D_mat2, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw =dlmread(filename,'\t',[43 1 50 nfiles]);
            elseif (strfind(D_mat3, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw =dlmread(filename,'\t',[92 1 99 nfiles]);
            elseif (strfind(D_mat4, 'Conc (SSFDPM):')) == 1
                broadband_fit_data_raw =dlmread(filename,'\t',[56 1 63 nfiles]);
            end
            
            broadband_fit_data=broadband_fit_data_raw([1:4 6:8],:);
            fclose(fid);
            old_concentrations = broadband_fit_data_raw(1:5,:);
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
                if data_fitting == 1
                   PatientData.(PointLabel).muaSpectrum(4,:)=data_fit(:,j);
                else
                   fdpmchrom3 = physioSetup(physio, wvs, 1);
                   muas = data(:,j);
                   [~, fdpmfit3(j).fitmua] = physioFit3(physio.fittype, fdpmchrom3, muas', muas', 0, old_concentrations(1:5,j)');
                   mua_fit_spec_old = fdpmfit3(j).fitmua;
                   PatientData.(PointLabel).muaSpectrum(4,:) = mua_fit_spec_old;
%                    plot(wvs,mua_fit_spec,wvs,mua_fit_spec_old,wvs,data(:,j))
%                    legend('new fit','old fit', 'mua')
                end
                PatientData.(PointLabel).pointlabel = headers{1,1}([false logical(header_index)]);
                PatientData.(PointLabel).tag = DataTags(:,j);
                PatientData.(PointLabel).chromo(1,:) = fdpmfit(1,j).phy(1,:);
                PatientData.(PointLabel).chromo(2,:) = broadband_fit_data(:,j)'; 
                PatientData.(PointLabel).chromo(3,:) = broadband_fit_data(:,j)'-fdpmfit(1,j).phy(1,:); 
                PatientData.(PointLabel).chromo(4,:) = ((broadband_fit_data(:,j)'-fdpmfit(1,j).phy(1,:))./broadband_fit_data(:,j)')*100; 
            end
            OutputFile = strcat(folders(k+2).name,'_waves',num2str(waves,'%02d'),'.mat');
            cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Breast Results')
            save(OutputFile,'PatientData')
        end
end

%}