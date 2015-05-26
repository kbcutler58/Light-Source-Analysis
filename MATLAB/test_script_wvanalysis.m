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
                
                fdpmchrom = physioSetup(physio, testwvs,1);
                [fdpmfit(j).phy, fdpmfit(j).fitmua, fdpmfit(j).ss] = physioFit(physio.fittype, fdpmchrom, muas', muas', 0);
%                 physio.chrom.maxs = [.1;.1;1;1;.1;.1;0;0];
%                 muas=data(indexes,j*2);
                fdpmchrom2 = physioSetup(physio, wvs, 1);
                concentrations = fdpmfit(j).phy(1,1:5);
                [~, fdpmfit2(j).fitmua] = physioFit3(physio.fittype, fdpmchrom2, muas', muas', 0, concentrations);
                mua_fit_spec = fdpmfit2(j).fitmua;

                PointLabel = strcat('Point',num2str(j,'%03d'));
                PatientData.(PointLabel).wavelengths = testwvs;
                PatientData.(PointLabel).muaSpectrum(1,:)=wvs';
                PatientData.(PointLabel).muaSpectrum(2,:)=mua_fit_spec;
                PatientData.(PointLabel).muaSpectrum(3,:)=data(:,j);
                PatientData.(PointLabel).pointlabel = headers{1,1}([false logical(header_index)]);
                PatientData.(PointLabel).tag = DataTags(:,j);
                PatientData.(PointLabel).chromo(1,:) = fdpmfit(1,j).phy(1,:);
                PatientData.(PointLabel).chromo(2,:) = broadband_fit_data(:,j)'; 
                PatientData.(PointLabel).chromo(3,:) = broadband_fit_data(:,j)'-fdpmfit(1,j).phy(1,:); 
                PatientData.(PointLabel).chromo(4,:) = ((broadband_fit_data(:,j)'-fdpmfit(1,j).phy(1,:))./broadband_fit_data(:,j)')*100; 