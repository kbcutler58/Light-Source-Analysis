 function [phy, fitmua] = physioFit(whichfit, chrom, mua, dmua, verbose, conc)
%   Computes concentration matrix given data and extinction matrix
%%%byh Recovers concentrations of selected physiological properties using
%%%extinction coefficients and recovered optical absorption.  

Eorig=chrom.E;
mua=mua';
[mua, excise_idx] = excise_vector(mua);  %must remove NaN's (can we test first?)
if(sum(excise_idx))
    chrom.E(excise_idx,:) = [];
    dmua(excise_idx) = [];
end

conc = conc'./chrom.mult;

Einv=inv(chrom.E.'*chrom.E)*chrom.E.';

%%%byh Standard use is case 3, constraining the fit to positive values
switch whichfit
    case 1	%simple weights
        %In this simple weighting scheme, I have divided the value of mua by the error
        %recovered by the Lev-Mar algorithm. This acts as a simple weight since the
        %least-squares algorithm suffers a greater penalty for larger values, which in
        %this case translates into an error value small comapered to the fitted value.
        mua_temp= mua ./ dmua;		%scale mua values
        Escale=zeros(size(chrom.E));
        for j=1:size(chrom.E,1)						%scale extinction coefficients
            Escale(j,:)=chrom.E(j,:)./dmua(j);
        end
        % Invert the matrix from M x = b, with b as mua's, x as concentrations, M as extinctions
        % into x = ((MT M)MT)^-1 b  since M is not necessarily square (if M is square, then x = M^-1 b)
        E_temp=inv(Escale.'*Escale)*Escale';
        conc=(E_temp*mua_temp);	%1/mm t
    case 2		%simple SVD
        [Ue,Se,Ve]=svd(chrom.E,0);
        for jj=1:length(Se)
            Se(jj,jj)=1./Se(jj,jj);
        end;
        conc=Ve*Se*Ue'*mua;
    case 0 %simple least squares fit
        conc=(Einv*mua);			% 1/mm  for mua
    case 3 %constrained, but only positive values
%         disp('Minimum & Maximum constrained least-squares algorithm engaged');
        %%%byh Concentrations recovered here
%         conc=lsqlin(chrom.E,mua,[],[],[],[],chrom.mins);
    case 4  %constrained LSQ
        disp('Minimum & Maximum constrained least-squares algorithm engaged');
        conc=lsqlin(chrom.E,mua,[],[],[],[],chrom.mins,chrom.maxs);
end

fitmua = Eorig*conc;
% residual = (mua-chrom.E*conc);
% chi = sqrt(sum(residual.^2));
% 
% conc_err = zeros(size(Einv,1),1);
% for j=1:size(Einv,1)		%just saying that error in conc is proportional to error in mua
%    conc_err(j)=sqrt(sum((Einv(j,:).*residual').^2));  %shot in the dark
% end
% 
% %do mult to get right units
% conc = conc.*chrom.mult;
% conc_err = conc_err.*chrom.mult;
% 
% %compute sat and thc?
% %%%byh Some values like THC, TOI, Sat are calculated from combinations of
% %%%the recovered chomophore concentrations
% if chrom.saton
%    hbO2=conc(chrom.hbO2_i);
%    hb=conc(chrom.hb_i);
%    hbO2_err=conc_err(chrom.hbO2_i);
%    hb_err =conc_err(chrom.hb_i);
%    thc = hbO2 + hb;
%    sat = hbO2/thc*100;
%    thc_err = sqrt(hbO2_err^2 + hb_err^2);
%    sat_err = sqrt((hbO2_err/thc)^2+(thc_err*hbO2/thc^2)^2)*100;
%    conc = [conc; thc; sat];	
%    conc_err = [conc_err; thc_err; sat_err];
% end
% if chrom.toion
%     fat=conc(chrom.fat_i);
%     fat_err=conc_err(chrom.fat_i);
%     h2o=conc(chrom.h2o_i);
%     h2o_err=conc_err(chrom.h2o_i);
%     toi=h2o*hb/fat;
%     toi_err=sqrt((hb_err/hb)^2+(fat_err/fat)^2+(h2o_err/h2o)^2);
%     conc=[conc; toi];
%     conc_err = [conc_err; toi_err];
% end
% 
% if verbose
%     for i = 1:length(chrom.names),disp(sprintf('[%s]\t= %3.2f +/- %2.3f %s', chrom.names{i}, conc(i), conc_err(i), chrom.units{i})); end        
% end
%     
phy = [];