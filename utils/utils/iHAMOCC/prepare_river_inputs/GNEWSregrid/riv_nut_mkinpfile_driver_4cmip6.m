% J. Schwinger 20.09.2024 (based on scripts from S. Gao 30.03.16)
% readin nutrient fluxes from GNEWS2 and other sources and 
% interpolate to 0.5x0.5 degree runoff grid or (if regrid2model=true)
% to the BLOM ocean grid (specified by resol_tag) and write into NetCDF files


clear all


global use_octave

use_octave=true;

if use_octave
    pkg load netcdf
    import_netcdf
end



datapath = './RiverineData/';
datapath_dic = './RiverineData/UHH/CarbonLateralExports2000';
%scenario = 'c70';
scenario = 'c00';  % c00 files have been used in NorESM so far
regrid2model = false;
vstr = '20250220';

%resol_tag = {'tnx0.25v4'};
%resol_tag = {'tnx0.5v1'};
resol_tag = {'tnx1v4'};
%resol_tag = {'tnx1v2'};
%resol_tag = {'gnx1v1'};
%resol_tag = {'tnx2v1'};

% read in data
    
    % Read lat/lon of river mouth locations for each basin in GNEWS2
    filename_basin = [datapath '/GNEWs2/' scenario 'NEWS_basins_red.csv'];
    data_basin  = importdata(filename_basin,';',1);
    data0 = data_basin.data(:,2:3);  % lon/lat of 6081 rivermouth

    % Read hydrology information of each river mouth in GNEWS2
    filename_hydro = [datapath '/GNEWs2/' scenario 'NEWS_hydrology.csv'];
    data_hydrol = importdata(filename_hydro,';',1);
    data1 = data_hydrol.data;
    
    % Read nutrient export data
    filename_export = [datapath '/GNEWs2/' scenario 'NEWS_river_exports.csv'];
    data_export = importdata(filename_export,';',1);
    data2 = data_export.data;

    % Read DIC export data        
    filename_dic = [datapath_dic '/NEWSOutputExports_2000Carbon_H2009.csv'];
    data_dic_r00 = importdata(filename_dic,',',1);
    dic_r00 = data_dic_r00.data(:,5);

    % Scale DIC with runoff: dic_r00/Qact_r00*Qact_{senario};
    % Total DIC changes among senarios
    filename_Qactr00 = [datapath '/GNEWs2/r00NEWS_hydrology.csv']; 
    data_r00 = importdata(filename_Qactr00,';',1);
    
    Qact_r00 = data_r00.data(:,4);
    Qact = data1(:,4);
    fac_dic = Qact./Qact_r00;
    idx = isnan(fac_dic) | isinf(fac_dic);
    fac_dic(idx) = 0.;
    data3_dic = dic_r00.*fac_dic;

    % sum(dic_r00)
    % sum(data3_dic)
    % sum(Qact_r00)
    % sum(Qact)

    % Scale DFe with runoff:
    % Fetotal*(sum(Qact)/sum(Qact_c70))*(Qact/sum(Qact))=Fetotal*Qact/sum(Qact_c70)
    Fe_total = 1.45e6;  % Mg/yr = 1.45 Tg/yr (Chester,1990) as input for yr1970, total DFe changes among senarios
    filename_Qactc70 = [datapath '/GNEWs2/c70NEWS_hydrology.csv']; 
    data_c70 = importdata(filename_Qactc70,';',1);
    Qact_c70 = data_c70.data(:,4);
    data4_fe = data1(:,4)/sum(Qact_c70)*Fe_total;
    %data4_fe = data1(:,4)/sum(data1(:,4))*Fe_total;  % scale Fe by Qact, no change in total flux 


    % convert units of runoff-data 
    cfac_rof = 1.0e12/(365.*24.*60.*60.); % runoff km3/yr-->kg/m2/s (mass=density*V;water density=1000kg/m3)
    data1(:,4:5) = data1(:,4:5)*cfac_rof;

    %  1.'BASINID' 2.'GN_lon' 3.'GN_lat' };
    data_loc = [ data1(:,1) data0 ];
    %  4.'Qact' 5.'Qnat' 6.'Ld_DIN' 7.'Ld_DIP' 8.'Ld_DON' 9.'Ld_DOP' 10.'Ld_DOC' 11.'Ld_DSi' 12.'Ld_PN' 13.'Ld_PP' 14.'Ld_POC' 15.'Ld_DIC' 16. 'Ld_Fe'
    data_riv = [ data1(:,4:5) data2(:,11:19) data3_dic data4_fe];

    [data_rivred] = riv_nut_calc_redfield(data_riv);

   
    if regrid2model
	
        disp('Regriding to model grid')
    
        % Regrid to model grid given by 'resol_tag'
        [idm,jdm,data_ocegrid] = regrid_riv_model(data_loc(:,2),data_loc(:,3),data_rivred,scenario,resol_tag);
    
        % Write data to file
        riv_nut_ncwrite_model(idm,jdm,scenario,data_ocegrid,resol_tag,vstr);

    else
	
        disp('Regriding to runoff grid')
	
        % Regrid to runoff grid
        [idm,jdm,lon,lat,data_rofgr] = regrid_riv_rofgr(data_loc(:,2),data_loc(:,3),data_rivred,scenario,resol_tag);
        
        % Write data to a file on the 0.5x0.5 runoff grid
        riv_nut_ncwrite_rofgr(idm,jdm,lon,lat,data_rofgr,scenario,vstr);

   end


