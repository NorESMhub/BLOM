% 07.04.16 sgao adapted from Riv_nutr_inp_wrtout.m
% calculate riverine input nutrients according to Redfield ratio and
% combine some nutrients

function [data_rivred] = riv_nut_calc_redfield(data_riv)

% unit convertion
facN = 1E6/14./1E3;
facP = 1E6/31./1E3;
facC = 1E6/12./1E3;
facSi = 1E6/28./1E3;
facFe = 1E6/55.8/1E3;   

% 1.'Qact' 2.'Qnat' 3.'DIN' 4.'DIP' 5.'DON' 6.'DOP' 7.'DOC' 8.'DSi' 9.'PN' 10.'PP' 11.'POC' 12.'DIC' 13.'Fe'};
fac = [1. 1. facN facP facN facP facC facSi facN facP facC facC facFe]; 

nx = size(data_riv,1);
riv_2d = zeros(nx,13);
for i=1:nx,
        riv_2d(i,:) = data_riv(i,:).*fac;
end

% calculate redfield ratio
rnit = 16.;
rcar = 122.;
riv_DON_cvt = riv_2d(:,5)/rnit;
riv_DOC_cvt = riv_2d(:,7)/rcar;
riv_PN_cvt = riv_2d(:,9)/rnit;
riv_POC_cvt = riv_2d(:,11)/rcar;
for i=1:nx,
    riv_idoc2d(i) = min([riv_2d(i,6),riv_DON_cvt(i),riv_DOC_cvt(i)]);
    res_DOP = max( riv_2d(i,6) - riv_idoc2d(i),0.);
    res_DON = max( riv_DON_cvt(i) - riv_idoc2d(i),0.);
    res_DOC = max( riv_DOC_cvt(i) - riv_idoc2d(i),0.);
    riv_idet2d(i) = min([riv_2d(i,10),riv_PN_cvt(i),riv_POC_cvt(i)]);
    res_PP = max( riv_2d(i,10) - riv_idet2d(i),0.);
    res_PN = max( riv_PN_cvt(i) - riv_idet2d(i),0.);
    res_POC = max( riv_POC_cvt(i) - riv_idet2d(i),0.);
    riv_DIP2d(i) = riv_2d(i,4) + res_DOP + res_PP;
    riv_DIN2d(i) = riv_2d(i,3) + res_DON*rnit + res_PN*rnit;
    riv_DIC2d(i) = riv_2d(i,12) + res_DOC*rcar + res_POC*rcar;
end

%              Qact Qnat      DIN        DIP        DSi         DIC        DET         DOC         DFe
data_rivred = [riv_2d(:,1:2), riv_DIN2d',riv_DIP2d',riv_2d(:,8),riv_DIC2d',riv_idet2d',riv_idoc2d',riv_2d(:,13)];


end
