clear all
close all
clc
rng('shuffle');
load(['../eqs_NZ_augmented_2025_fullhik_nolocuncert.mat']);
tic
wgs84 = wgs84Ellipsoid("km");

%%
maxdepth = 40;
loncens = [floor(min(eqs_AC.poly50_lons)*5)/5:0.1:ceil(max(eqs_AC.poly50_lons)*5)/5];
loncens = loncens(1:end-1)/2 + loncens(2:end)/2;
latcens = [floor(min(eqs_AC.poly50_lats)*5)/5:0.1:ceil(max(eqs_AC.poly50_lats)*5)/5];
latcens = latcens(1:end-1)/2 + latcens(2:end)/2;
[longrid,latgrid] = meshgrid(loncens,latcens);

depthspacing_1km = 1;
depthbinedges_1km = [0:depthspacing_1km:maxdepth]';
depthcens_1km = depthbinedges_1km(1:end-1) + depthspacing_1km/2;
[~,~,depthgrid_1km] = meshgrid(loncens,latcens,depthcens_1km);
density_3D = NaN*depthgrid_1km;
%return

in50 = inpolygon(longrid,latgrid,eqs_AC.poly50_lons,eqs_AC.poly50_lats);

load(['../hik_kerm_fault_300km_wgs84_poslon.txt'])
hik_CW.lons = hik_kerm_fault_300km_wgs84_poslon(:,1);
hik_CW.lats = hik_kerm_fault_300km_wgs84_poslon(:,2);
hik_CW.depths = -hik_kerm_fault_300km_wgs84_poslon(:,3)/1e3;

hik_CW.interp = scatteredInterpolant(hik_CW.lons,hik_CW.lats,hik_CW.depths,'linear','none');
hikdepthgrid = hik_CW.interp(longrid,latgrid);

hik_CW.boundary = boundary(hik_CW.lons,hik_CW.lats,0.75);
hik_outerrise_lons = [hik_CW.lons(hik_CW.boundary(18:150)); 195; hik_CW.lons(hik_CW.boundary(18)) + 40*(hik_CW.lons(hik_CW.boundary(18)) - hik_CW.lons(hik_CW.boundary(17))); hik_CW.lons(hik_CW.boundary(18))];
hik_outerrise_lats = [hik_CW.lats(hik_CW.boundary(18:150)); hik_CW.lats(hik_CW.boundary(150)); hik_CW.lats(hik_CW.boundary(18)) + 50*(hik_CW.lats(hik_CW.boundary(18)) - hik_CW.lats(hik_CW.boundary(17)));; hik_CW.lats(hik_CW.boundary(18))];

in50(inpolygon(longrid,latgrid,hik_outerrise_lons,hik_outerrise_lats)) = 0;

puy_HS.interp = scatteredInterpolant(puy_HS.lons,puy_HS.lats,puy_HS.depths,'linear','none');
puydepthgrid = puy_HS.interp(longrid,latgrid);

puy_HS.boundary = boundary(puy_HS.lons,puy_HS.lats,1);
puy_outerrise_lons = [puy_HS.lons(puy_HS.boundary(553:582)); puy_HS.lons(puy_HS.boundary(1:276)); puy_HS.lons(puy_HS.boundary(276)) + 2*sind(225); 162; puy_HS.lons(puy_HS.boundary(553)) + 2*sind(315); puy_HS.lons(puy_HS.boundary(553))];
puy_outerrise_lats = [puy_HS.lats(puy_HS.boundary(553:582)); puy_HS.lats(puy_HS.boundary(1:276)); puy_HS.lats(puy_HS.boundary(276)) + 2*cosd(225); -40; puy_HS.lats(puy_HS.boundary(553)) + 2*cosd(315); puy_HS.lats(puy_HS.boundary(553))];
in50(inpolygon(longrid,latgrid,puy_outerrise_lons,puy_outerrise_lats)) = 0;
clear eqs_AC
sum(in50(:))

%%
depthmodel = readtable(['../upperplate_depthdensitypercell_v0pt95_locuncert10_01-Nov-2025.dat']);
depthmodel.lons = depthmodel.Var1;
depthmodel.lats = depthmodel.Var2;
depthmodel.depthcens = depthmodel.Var3;
depthmodel.density = depthmodel.Var4;

ratemodel = readtable(['upperplate_2Dratemap_v0pt9_gaulap_lastcolumnsup_wmed_24-Oct-2025.dat']);
ratemodel.lons = ratemodel.Var1;
ratemodel.lats = ratemodel.Var2;
ratemodel.density = ratemodel.Var3;

%%
for lonindex = 1:length(loncens)
     lonindex
     ratemodel_lonindex = abs(ratemodel.lons - loncens(lonindex))<1e-3;
     depthmodel_lonindex = abs(depthmodel.lons - loncens(lonindex))<1e-3;
     for latindex = 1:length(latcens)
          if in50(latindex,lonindex)
          ratemodelindex = and(ratemodel_lonindex,abs(ratemodel.lats - latcens(latindex))<1e-3);
          depthmodelindex_2D = and(depthmodel_lonindex,abs(depthmodel.lats - latcens(latindex))<1e-3);
          for depthindex = 1:length(depthcens_1km)
               depthmodelindex = and(depthmodelindex_2D,abs(depthmodel.depthcens - depthcens_1km(depthindex))<1e-3);
               density_3D(latindex,lonindex,depthindex) = ratemodel.density(ratemodelindex)*depthmodel.density(depthmodelindex);
          end
          end
     end
end

%%
load colormaps.mat
coastlines = readtable(['nz-coastline-mean-high-water_180.dat']);
lakes = readtable(['nz-lake-polygons-topo-1500k.dat']);
load(['cfm_1pt0_vertices_NaNs.txt']);
cfm_1pt0_vertices_NaNs(cfm_1pt0_vertices_NaNs(:,1)<0,1) = cfm_1pt0_vertices_NaNs(cfm_1pt0_vertices_NaNs(:,1)<0,1) + 360;
%%
clim_prc = 5e-5;
%depthcens = depthcens;
for depthindex = 1:length(depthcens_1km)
     %     plotdepth = depthcens(plotindex);
     %     densitygrid_toplot = densitygrid_wmed(:,:,plotindex);
     figure(depthindex); clf; hold on; box on;
     set(gca,'DataAspectRatio',[1/cosd(-41) 1 111.1/5])
     set(gcf,'position',[0 725 525 725])
     %plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.8 0.8 0.8],'linewidth',0.25)
     %plot(eqs_AC.poly50_lons,eqs_AC.poly50_lats,'color','k','linewidth',0.5)
     plot(coastlines.X,coastlines.Y,'color',[0 0.25 0.375],'linewidth',0.5)
     plot(lakes.X,lakes.Y,'color',[0 0.25 0.375],'linewidth',0.5)

     scatter3(hik_CW.lons,hik_CW.lats,-hik_CW.depths,12.5,hik_CW.depths,'v','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
     scatter3(puy_HS.lons,puy_HS.lats,-puy_HS.depths,2.5,puy_HS.depths,'^','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
     %plot(hik_CW.lons,hik_CW.lats,'*')
     %text(hik_CW.lons,hik_CW.lats,num2str([1:length(hik_CW.lats)]'))
     pcolor(longrid,latgrid,log10(density_3D(:,:,depthindex))); shading interp;

     xlim([165.5 179.5])
     ylim([-48 -34])
     caxis([-6 -4])
%     caxis([0 0.15])
     zlim([-40 0])
     colormap((jet4_1000))
     hcb = colorbar;
     hcb.Position = [0.07 0.4675 0.04 0.45];
     hcb.AxisLocation = 'in';
     %hcb.Direction = "reverse";
     fontname(gcf,'Minion Pro')
     %fontweight(gcf,'bold')
     %set(gcf,'position',[0 800 800 1010])
     ax = gca;
     ax.FontSize = 13.5;
     outerpos = ax.OuterPosition;
     ti = ax.TightInset;
     left = outerpos(1) + ti(1);
     bottom = outerpos(2) - 0.02;
     ax_width = outerpos(3) - ti(1)*1.2;
     ax_height = outerpos(4) - ti(2) - ti(4);
     ax.Position = [left bottom ax_width ax_height];
     %return
     %return
     title({['3D density of upper-plate earthquakes at ' num2str(depthcens_1km(depthindex)) ' km depth'],['(2D rate model x 3D depth density model)']},'fontsize',20,'fontweight','normal')
     %return
     print(depthindex,['ratemodel3D_upperplate_' num2str(depthcens_1km(depthindex)) '_' date '_log10.png'],'-dpng','-r280');
end