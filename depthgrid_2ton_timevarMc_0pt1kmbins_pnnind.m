clear all
close all
clc
rng('shuffle');
addpath(['~/My Drive (john.c.rollins@gmail.com)//utils/'])
load(['~/My Drive (john.c.rollins@gmail.com)//NZ/Material/seismicity/eqs_NZ_augmented_2025_fullhik_nolocuncert.mat']);
tic
wgs84 = wgs84Ellipsoid("km");

%%
numperiods = 500;
maxdepth = 40;
startdate = 1964;
enddate = 2021;
minbufferlength = 5;
maxbufferlength = 10;
minnum = 32;
minduration = 1;
maxMc = 4.5;
locuncert = 10;
%minmag = 2.8;

%minmag = -1;
sel.ind = or(eqs_AC.depths_uncorr~=eqs_AC.depths,strcmp(eqs_AC.method,'NonLinLoc'));
sel.ind = and(sel.ind,eqs_AC.depths_fixed_anduncorrected==0);
sel.ind = and(sel.ind,and(eqs_AC.depthuncert<14.999,eqs_AC.depthuncert>=0.05));
sel.ind(and(isnan(eqs_AC.depths_uncorr),round(eqs_AC.depths/5)==eqs_AC.depths/5)) = 0;
sel.ind = and(sel.ind,eqs_AC.crustal_weights==1);
sel.ind = and(sel.ind,and(and(eqs_AC.dates>=startdate,eqs_AC.dates<enddate),eqs_AC.depths<=maxdepth));
sel.ind = and(sel.ind,eqs_AC.poly100filter);
%sel.ind = and(sel.ind,eqs_AC.Mpref>=minmag);
%sel.ind = and(sel.ind,inpolygon(eqs_AC.lons,eqs_AC.lats,bigpolygon.lons,bigpolygon.lats));

%sel.ind = and(sel.ind,eqs_AC.Mpref>=minmag);
normconst = 1./(sqrt(2.*pi).*eqs_AC.depthuncert);
%labelmag = 6.95;

sel.depths = eqs_AC.depths(sel.ind);
sel.lons = eqs_AC.lons(sel.ind);
sel.lats = eqs_AC.lats(sel.ind);
sel.depthuncert = eqs_AC.depthuncert(sel.ind);
sel.normconst = normconst(sel.ind);
sel.dates = eqs_AC.dates(sel.ind);
sel.in50 = eqs_AC.poly50filter(sel.ind);
%sel.in100 = eqs_AC.poly100filter(sel.ind);
sel.mags = eqs_AC.Mpref(sel.ind);
sel.ind = [];
minmag = min(sel.mags);

%return
%sel.weights = eqs_AC.weights(sel.ind);

%return
%%
loncens = [floor(min(eqs_AC.poly50_lons)*5)/5:0.1:ceil(max(eqs_AC.poly50_lons)*5)/5];
loncens = loncens(1:end-1)/2 + loncens(2:end)/2;
latcens = [floor(min(eqs_AC.poly50_lats)*5)/5:0.1:ceil(max(eqs_AC.poly50_lats)*5)/5];
latcens = latcens(1:end-1)/2 + latcens(2:end)/2;
[longrid,latgrid] = meshgrid(loncens,latcens);

depthspacing_fine = 0.025;
depthbinedges_fine = [0:depthspacing_fine:maxdepth]';
depthcens_fine = depthbinedges_fine(1:end-1) + depthspacing_fine/2;

depthspacing_1km = 1;
depthbinedges_1km = [0:depthspacing_1km:maxdepth]';
depthcens_1km = depthbinedges_1km(1:end-1) + depthspacing_1km/2;
[~,~,depthgrid_1km] = meshgrid(loncens,latcens,depthcens_1km);
%return

sel.pdfs_fine = sel.normconst.*exp((-(sel.depths - depthcens_fine').^2)./(2.*sel.depthuncert.^2));
in50 = inpolygon(longrid,latgrid,eqs_AC.poly50_lons,eqs_AC.poly50_lats);

load(['~/My Drive (john.c.rollins@gmail.com)//NZ/Material/Hikurangi/hik_kerm_fault_300km_wgs84_poslon.txt'])
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
%return
%%
deepestbincenind_fine = NaN*longrid;
maxdepth_fine = NaN*longrid;
deepestbincenind_1km = NaN*longrid;
maxdepth_1km = NaN*longrid;
inthisdepthbin_1km = cell(length(latcens),length(loncens));

for lonindex = 1:length(loncens)
     for latindex = 1:length(latcens)
          if in50(latindex,lonindex)

               insubd_fine = or(depthbinedges_fine>hikdepthgrid(latindex,lonindex),depthbinedges_fine>puydepthgrid(latindex,lonindex));
%               incrust_fine = insubd_fine==0;
               if sum(insubd_fine)>0
                    deepestbincenind_fine(latindex,lonindex) = find(insubd_fine==0,1,'last')-1;
               else
                    deepestbincenind_fine(latindex,lonindex) = length(depthcens_fine);
               end
               maxdepth_fine(latindex,lonindex) = depthbinedges_fine(deepestbincenind_fine(latindex,lonindex) + 1);

               deepestbincenind_1km(latindex,lonindex) = find(depthcens_1km<ceil(maxdepth_fine(latindex,lonindex)),1,'last');

               for depthind = 1:deepestbincenind_1km(latindex,lonindex)
                    inthisdepthbin_1km{latindex,lonindex}{depthind} = find(and(and(...
                         depthcens_fine>depthbinedges_1km(depthind),...
                         depthcens_fine<depthbinedges_1km(depthind+1)),depthcens_fine<=maxdepth_fine(latindex,lonindex)));
               end

               if maxdepth_fine(latindex,lonindex) == 0
                    in50(latindex,lonindex) = 0;
               end
          end
     end
end
sum(in50(:))

%%
train = [];
test = [];
nnfracarray = 2.^([-15:0]');
nnsize = length(nnfracarray);
densityfun_1km_out = cell(numperiods,1);

logprob_fine = zeros(numperiods,nnsize);
logprob_sup_fine = zeros(numperiods,nnsize);

aicc_test = zeros(numperiods,nnsize);
aicc_test_sup = zeros(numperiods,nnsize);
worsethansup = ones(numperiods,nnsize)==0;

sumincell_test_fine_all = zeros(numperiods,nnsize);

Mc = NaN(numperiods,1);
bufferlength = Mc;
%densityvecs_all = cell(numperiods,1);
for p = 1:numperiods
     p
     Mc(p) = (minmag - 1e-3) + rand(1)*(maxMc - minmag);
     if Mc(p)>=3.8
          startdate_local = 1964;
     elseif Mc(p)>=3.45
          startdate_local = 1987;
     else
          startdate_local = 2000;
     end
     bufferlength(p) = minbufferlength + rand(1)*(maxbufferlength - minbufferlength);
     issue = 1;

     while issue==1
          issue = 0;

          if rand(1)>=0.5
               divdate = rand(1)*(enddate - startdate_local) + startdate_local;

               if rand(1)>=0.5
                    train.startdate1(p) = startdate_local;
                    train.enddate1(p) = divdate - bufferlength(p)/2;
                    train.startdate2(p) = NaN;
                    train.enddate2(p) = NaN;

                    test.startdate1(p) = divdate + bufferlength(p)/2;
                    test.enddate1(p) = enddate;
                    test.startdate2(p) = NaN;
                    test.enddate2(p) = NaN;
               else
                    train.startdate1(p) = divdate + bufferlength(p)/2;
                    train.enddate1(p) = enddate;
                    train.startdate2(p) = NaN;
                    train.enddate2(p) = NaN;

                    test.startdate1(p) = startdate_local;
                    test.enddate1(p) = divdate - bufferlength(p)/2;
                    test.startdate2(p) = NaN;
                    test.enddate2(p) = NaN;
               end

               train.duration(p) = train.enddate1(p) - train.startdate1(p);
               train.ind{p} = and(and(sel.dates>=train.startdate1(p),sel.dates<=train.enddate1(p)),sel.mags>=Mc(p));
               test.duration(p) = test.enddate1(p) - test.startdate1(p);
               test.ind{p} = and(and(sel.dates>=test.startdate1(p),sel.dates<=test.enddate1(p)),sel.mags>=Mc(p));

          else
               divdates = rand(2,1)*(enddate - startdate_local) + startdate_local;

               if rand(1)>=0.5
                    train.startdate1(p) = divdates(1) + bufferlength(p)/2;
                    train.enddate1(p) = divdates(2) - bufferlength(p)/2;
                    train.startdate2(p) = NaN;
                    train.enddate2(p) = NaN;

                    train.duration(p) = train.enddate1(p) - train.startdate1(p);
                    train.ind{p} = and(and(sel.dates>=train.startdate1(p),sel.dates<=train.enddate1(p)),sel.mags>=Mc(p));

                    test.startdate1(p) = startdate_local;
                    test.enddate1(p) = divdates(1) - bufferlength(p)/2;
                    test.startdate2(p) = divdates(2) + bufferlength(p)/2;
                    test.enddate2(p) = enddate;

                    test.duration(p) = (test.enddate1(p) - test.startdate1(p)) + (test.enddate2(p) - test.startdate2(p));
                    test.ind{p} = and(or(and(sel.dates>=test.startdate1(p),sel.dates<=test.enddate1(p)),and(sel.dates>=test.startdate2(p),sel.dates<=test.enddate2(p))),sel.mags>=Mc(p));
               else
                    test.startdate1(p) = divdates(1) + bufferlength(p)/2;
                    test.enddate1(p) = divdates(2) - bufferlength(p)/2;
                    test.startdate2(p) = NaN;
                    test.enddate2(p) = NaN;

                    test.duration(p) = test.enddate1(p) - test.startdate1(p);
                    test.ind{p} = and(and(sel.dates>=test.startdate1(p),sel.dates<=test.enddate1(p)),sel.mags>=Mc(p));

                    train.startdate1(p) = startdate_local;
                    train.enddate1(p) = divdates(1) - bufferlength(p)/2;
                    train.startdate2(p) = divdates(2) + bufferlength(p)/2;
                    train.enddate2(p) = enddate;

                    train.duration(p) = (train.enddate1(p) - train.startdate1(p)) + (train.enddate2(p) - train.startdate2(p));
                    train.ind{p} = and(or(and(sel.dates>=train.startdate1(p),sel.dates<=train.enddate1(p)),and(sel.dates>=train.startdate2(p),sel.dates<=train.enddate2(p))),sel.mags>=Mc(p));
               end
          end

          if or(abs(train.startdate1 - test.enddate1)<bufferlength(p),abs(train.startdate2 - test.enddate2)<bufferlength(p))
               issue = 1;
          end
          if or(abs(train.startdate2 - test.enddate1)<bufferlength(p),abs(train.startdate1 - test.enddate2)<bufferlength(p))
               issue = 1;
          end
          if sum(and(train.ind{p},sel.in50))<minnum
               issue = 1;
          end
          if sum(and(test.ind{p},sel.in50))<minnum
               issue = 1;
          end
          if or(test.duration(p)<minduration,train.duration(p)<minduration)
               issue = 1;
          end
     end
     train.weightsum_50(p) = sum(and(train.ind{p},sel.in50));

end
%return
%%
parfor p = 1:numperiods
     trainlons = sel.lons(train.ind{p}) + randn(sum(train.ind{p}),1).*locuncert./(111.1.*cosd(sel.lats(train.ind{p})));
     trainlats = sel.lats(train.ind{p}) + randn(sum(train.ind{p}),1)*locuncert/111.1;
     trainpdfs_fine = sel.pdfs_fine(train.ind{p},:);
     trainweightsum_50 = sum(and(train.ind{p},sel.in50));
     numnearest_local = round(nnfracarray*trainweightsum_50);

     testlons = sel.lons(test.ind{p}) + randn(sum(test.ind{p}),1).*locuncert./(111.1.*cosd(sel.lats(test.ind{p})));
     testlats = sel.lats(test.ind{p}) + randn(sum(test.ind{p}),1)*locuncert/111.1;
     testdepths = sel.depths(test.ind{p}) + randn(sum(test.ind{p}),1).*sel.depthuncert(test.ind{p});

     for nnind = 1:nnsize
          densityfun_1km_out{p,nnind} = cell(size(longrid,1),size(longrid,2));
     end

     for lonindex = 1:length(loncens)
          for latindex = 1:length(latcens)
               if in50(latindex,lonindex)
                    [p lonindex latindex]
                    thislat = latcens(latindex);
                    thislon = loncens(lonindex);

                    inthiscell_fine = and(and(...
                         and(testlons>=thislon-0.05,testlons<thislon+0.05),...
                         and(testlats>=thislat-0.05,testlats<thislat+0.05)),testdepths<=maxdepth_fine(latindex,lonindex));
                    suminthiscell_fine = sum(inthiscell_fine);
                    inthisdepthbin_1km_thisll = inthisdepthbin_1km{latindex,lonindex};

                    % inthiscell_1km = and(and(...
                    %      and(testlons>=thislon-0.05,testlons<thislon+0.05),...
                    %      and(testlats>=thislat-0.05,testlats<thislat+0.05)),testdepths<=maxdepth_1km(latindex,lonindex));
                    % suminthiscell_1km = sum(inthiscell_1km);

                    % inthiscell_2km = and(and(...
                    %      and(testlons>=thislon-0.05,testlons<thislon+0.05),...
                    %      and(testlats>=thislat-0.05,testlats<thislat+0.05)),testdepths<=maxdepth_2km(latindex,lonindex));
                    % suminthiscell_2km = sum(inthiscell_2km);
                     
                    deepestbincenind_fine_thisll = deepestbincenind_fine(latindex,lonindex);
                    deepestbincenind_1km_thisll = deepestbincenind_1km(latindex,lonindex);
                    maxdepth_fine_thisll = maxdepth_fine(latindex,lonindex);
                    %inthisdepthbin_1km = cell(deepestbincenind_1km_thisll,1);

                    if suminthiscell_fine>0
                         suponitsside_fine = ones(deepestbincenind_fine_thisll,1)/deepestbincenind_fine_thisll;

                         numperdepth_fine = histcounts(testdepths(inthiscell_fine),depthbinedges_fine(1:deepestbincenind_fine_thisll+1))';
                         prednum_sup_fine = suponitsside_fine*suminthiscell_fine;

                         factterm_fine = log(factorial(numperdepth_fine));
                         badinds_fine = find(factterm_fine==Inf);
                         for badind = 1:length(badinds_fine)
                              factterm_fine(badinds_fine(badind)) = sum(log([1:numperdepth_fine(badinds_fine(badind))]));
                         end

                         logprob_sup_fine_thisll = -prednum_sup_fine...
                              + numperdepth_fine.*log(prednum_sup_fine)...
                              - factterm_fine;
                    end

                    % deepestbincenind_1km_thisll = deepestbincenind_1km(latindex,lonindex);
                    % if suminthiscell_1km>0
                    %      suponitsside_1km = ones(deepestbincenind_1km_thisll,1)/deepestbincenind_1km_thisll;
                    % 
                    %      numperdepth_1km = histcounts(testdepths(inthiscell_1km),depthbinedges_1km(1:deepestbincenind_1km_thisll+1))';
                    %      prednum_sup_1km = suponitsside_1km*suminthiscell_1km;
                    % 
                    %      factterm_1km = log(factorial(numperdepth_1km));
                    %      badinds_1km = find(factterm_1km==Inf);
                    %      for badind = 1:length(badinds_1km)
                    %           factterm_1km(badinds_1km(badind)) = sum(log([1:numperdepth_1km(badinds_1km(badind))]));
                    %      end
                    % 
                    %      logprob_sup_1km_thisll = -prednum_sup_1km...
                    %           + numperdepth_1km.*log(prednum_sup_1km)...
                    %           - factterm_1km;
                    % end

                    eqcelldists = distance(trainlats,trainlons,thislat,thislon,wgs84);
                    [~,sortorder] = sort(eqcelldists);

                    densityfun_prev_fine = zeros(deepestbincenind_fine_thisll,1);


%                    densityfun_prev_1km = zeros(deepestbincenind_1km_thisll,1);
%                    densityfun_prev_2km = zeros(deepestbincenind_2km_thisll,1);

                    for nnind = 1:nnsize
                         if nnind>1
                              nearby_new = sortorder((numnearest_local(nnind-1)+1):numnearest_local(nnind));
                         else
                              nearby_new = sortorder(1:numnearest_local(nnind));
                         end

                         densityfun_fine = densityfun_prev_fine + sum(trainpdfs_fine(nearby_new,[1:deepestbincenind_fine_thisll]),1)';
                         densityfun_prev_fine = densityfun_fine;
                         if sum(densityfun_fine)>0
                              densityfun_fine_norm = densityfun_fine/sum(densityfun_fine);
                         else
                              densityfun_fine_norm = densityfun_fine;
                         end

                         if suminthiscell_fine>0
                              prednum_fine = densityfun_fine_norm*suminthiscell_fine;
                              logprob_fine_thisll = -prednum_fine...
                                   + numperdepth_fine.*log(prednum_fine)...
                                   - factterm_fine;
                              logprob_fine_thisll(and(prednum_fine==0,numperdepth_fine==0)) = 0;
                              logprob_fine(p,nnind) = logprob_fine(p,nnind) + sum(logprob_fine_thisll);
                              logprob_sup_fine(p,nnind) = logprob_sup_fine(p,nnind) + sum(logprob_sup_fine_thisll);
                              sumincell_test_fine_all(p,nnind) = sumincell_test_fine_all(p,nnind) + suminthiscell_fine;
                         end

                         % densityfun_1km = densityfun_prev_1km + sum(trainpdfs_1km(nearby_new,[1:deepestbincenind_1km_thisll]),1)';
                         % densityfun_prev_1km = densityfun_1km;
                         % if sum(densityfun_1km)>0
                         %      densityfun_1km_norm = densityfun_1km/sum(densityfun_1km);
                         % else
                         %      densityfun_1km_norm = densityfun_1km;
                         % end
                         for depthind = 1:deepestbincenind_1km_thisll
                              densityfun_1km_out{p,nnind}{latindex,lonindex}(depthind) = sum(densityfun_fine_norm(inthisdepthbin_1km_thisll{depthind}));
                         end
                         %return

                         % if suminthiscell_1km>0
                         %      prednum_1km = densityfun_1km_norm*suminthiscell_1km;
                         %      logprob_1km_thisll = -prednum_1km...
                         %           + numperdepth_1km.*log(prednum_1km)...
                         %           - factterm_1km;
                         %      logprob_1km_thisll(and(prednum_1km==0,numperdepth_1km==0)) = 0;
                         %      logprob_1km(p,nnind) = logprob_1km(p,nnind) + sum(logprob_1km_thisll);
                         %      logprob_sup_1km(p,nnind) = logprob_sup_1km(p,nnind) + sum(logprob_sup_1km_thisll);
                         %      sumincell_test_1km_all(p,nnind) = sumincell_test_1km_all(p,nnind) + suminthiscell_1km;
                         % end


                    end
               end
          end
     end

     for nnind = 1:nnsize
         aicc_test(p,nnind) = -2*logprob_fine(p,nnind) + 2 + 2./(sumincell_test_fine_all(p,nnind) - 2);
         aicc_test_sup(p,nnind) = -2*logprob_sup_fine(p,nnind) + 2 + 2./(sumincell_test_fine_all(p,nnind) - 2);

         if aicc_test(p,nnind)>aicc_test_sup(p,nnind)
              worsethansup(p,nnind) = 1;
         end
     end
end

delete(gcp('nocreate'))

save(['upperplate_depthmodel_' num2str(startdate) 'to' num2str(enddate) '_md' num2str(maxdepth) '_timevarMc_min' num2str(minmag) '_' num2str(numperiods) 'periods_minnum' num2str(minnum) '_0pt25kmbins_' date '.mat'],'-v7.3')
%%
badones = or(aicc_test==Inf,isnan(aicc_test));
worsethansup(badones) = 1;
aicc_test(badones) = Inf;
allbad = sum(worsethansup==1,2)==nnsize;
sum(allbad)
aicc_test(allbad,:) = Inf;

ilpe = aicc_test./(2*sumincell_test_fine_all);
ilpe(isnan(ilpe)) = Inf;
prob_model = exp(min(ilpe(:)) - ilpe);
prob_hyper = sum(prob_model);

prob_final = repmat(prob_hyper,numperiods,1);
prob_final(allbad) = 0;
prob_final = prob_final/sum(prob_final(:));
prob_final_vec = prob_final(:);

%%
depthcens_length = length(depthcens_1km);
depthdensitypercell_wmed = NaN*depthgrid_1km;
depthdensitypercell_w2 = NaN*depthgrid_1km;
depthdensitypercell_w16 = NaN*depthgrid_1km;
depthdensitypercell_w84 = NaN*depthgrid_1km;
depthdensitypercell_w98 = NaN*depthgrid_1km;
densityfun_background = NaN(depthcens_length,1);

%%
for lonindex = 1:length(loncens)
     %     lonindex
     for latindex = 1:length(latcens)
          if in50(latindex,lonindex)
               deepestbincenind_1km_thisll = deepestbincenind_1km(latindex,lonindex);
               densityvecs_thisll = NaN(numperiods,nnsize,depthcens_length);
               [lonindex latindex]
               for p = 1:numperiods
                    %                    if nnz(worsethansup(p,:)>0)
                    for nnind = 1:nnsize
                         if prob_final(p,nnind)>0
                              densityfuns_local = densityfun_1km_out{p,nnind}{latindex,lonindex};
                              densityfun_out = densityfun_background;
                              densityfun_out(1:deepestbincenind_1km_thisll) = densityfuns_local;
                              densityvecs_thisll(p,nnind,:) = densityfun_out;
                         end
                         %                         end
                    end
                    %                    end
               end

               %toc

               for depthindex = 1:depthcens_length
                    densityvecs_thisdepth = densityvecs_thisll(:,:,depthindex);
                    densityvecs_thisdepth = densityvecs_thisdepth(:);
                    if sum(~isnan(densityvecs_thisdepth))>0
                         [~,sortorder] = sort(densityvecs_thisdepth);

                         densityvecs_thisdepth = densityvecs_thisdepth(sortorder);
                         scores_sorted = prob_final_vec(sortorder);
                         scores_cumsum = cumsum(scores_sorted);
                         scores_sum = sum(scores_sorted);

                         depthdensitypercell_wmed(latindex,lonindex,depthindex) = densityvecs_thisdepth(find(scores_cumsum>=scores_sum/2,1,'first'));
                         depthdensitypercell_w16(latindex,lonindex,depthindex) = densityvecs_thisdepth(find(scores_cumsum>=scores_sum*0.15865,1,'first'));
                         depthdensitypercell_w84(latindex,lonindex,depthindex) = densityvecs_thisdepth(find(scores_cumsum>=scores_sum*0.84135,1,'first'));
                         depthdensitypercell_w2(latindex,lonindex,depthindex) = densityvecs_thisdepth(find(scores_cumsum>=scores_sum*0.02275,1,'first'));
                         depthdensitypercell_w98(latindex,lonindex,depthindex) = densityvecs_thisdepth(find(scores_cumsum>=scores_sum*0.97725,1,'first'));
                    end
               end
               %toc

          end
     end
end

%delete(gcp('nocreate'))

%%
d90 = NaN*longrid;
for lonindex = 1:length(loncens)
     lonindex
     for latindex = 1:length(latcens)
          depthvec_wmed = depthdensitypercell_wmed(latindex,lonindex,:);
          depthvec_sum = sum(depthvec_wmed(~isnan(depthvec_wmed)));
          depthdensitypercell_wmed(latindex,lonindex,:) = depthvec_wmed/depthvec_sum;
          depthdensitypercell_w2(latindex,lonindex,:) = depthdensitypercell_w2(latindex,lonindex,:)/depthvec_sum;
          depthdensitypercell_w16(latindex,lonindex,:) = depthdensitypercell_w16(latindex,lonindex,:)/depthvec_sum;
          depthdensitypercell_w84(latindex,lonindex,:) = depthdensitypercell_w84(latindex,lonindex,:)/depthvec_sum;
          depthdensitypercell_w98(latindex,lonindex,:) = depthdensitypercell_w98(latindex,lonindex,:)/depthvec_sum;

          if in50(latindex,lonindex)
               depthvec_wmed = depthdensitypercell_wmed(latindex,lonindex,:);
               depthvec_sum = sum(depthvec_wmed(~isnan(depthvec_wmed)));
               d90_withoutsubd = depthcens_1km(find(cumsum(depthvec_wmed)>=depthvec_sum*0.9,1,'first'));
               d90(latindex,lonindex) = min([d90_withoutsubd depthcens_1km(deepestbincenind_1km(latindex,lonindex))]);
          end
     end
end

%%
ffile = fopen(['upperplate_depthdensitypercell_v0pt95_locuncert' num2str(locuncert) '_' date '.dat'],'wt');
for lonindex = 1:length(loncens)
     lonindex
     for latindex = 1:length(latcens)
          for depthindex = 1:length(depthcens_1km)
               fprintf(ffile,'%.2f %.2f %.2f %.3e\n',loncens(lonindex),latcens(latindex),depthcens_1km(depthindex),depthdensitypercell_wmed(latindex,lonindex,depthindex));
          end
     end
end
fclose(ffile);

%%
ffile = fopen(['upperplate_d90_v0pt95_' date '.dat'],'wt');
for lonindex = 1:length(loncens)
     lonindex
     for latindex = 1:length(latcens)
          fprintf(ffile,'%.2f %.2f %.2f %.3e\n',loncens(lonindex),latcens(latindex),d90(latindex,lonindex));
     end
end
fclose(ffile);

%return
%%
%save(['upperplate_depthmodel_' num2str(startdate) 'to' num2str(enddate) '_md' num2str(maxdepth) '_allbins_timevarMc_' num2str(numperiods) 'periods_minnum' num2str(minnum) '_lu' num2str(locuncert) '_' date '.mat'],'-v7.3')

%%
load colormaps.mat
coastlines = readtable(['~/My Drive (john.c.rollins@gmail.com)//utils/NZcoastlines/nz-coastline-mean-high-water/nz-coastline-mean-high-water_180.dat']);
lakes = readtable(['~/My Drive (john.c.rollins@gmail.com)//utils/NZcoastlines/nz-lake-polygons-topo-1500k/nz-lake-polygons-topo-1500k.dat']);
load(['~/My Drive (john.c.rollins@gmail.com)/NZ/Material/faults/CFM_2022/cfm_1pt0_vertices_NaNs.txt']);
cfm_1pt0_vertices_NaNs(cfm_1pt0_vertices_NaNs(:,1)<0,1) = cfm_1pt0_vertices_NaNs(cfm_1pt0_vertices_NaNs(:,1)<0,1) + 360;
clim_prc = prctile(depthdensitypercell_wmed(:),99.865);
%%
%depthcens = depthcens;
for plotindex = 1:length(depthcens_1km)
     %     plotdepth = depthcens(plotindex);
     %     densitygrid_toplot = densitygrid_wmed(:,:,plotindex);
     figure(plotindex); clf; hold on; box on;
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
     pcolor(longrid,latgrid,depthdensitypercell_wmed(:,:,plotindex)); shading interp;

     xlim([165.5 179.5])
     ylim([-48 -34])
%     caxis([0 clim_prc])
     caxis([0 0.15])
     zlim([-40 0])
     colormap((jet4w_1000))
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
     title({['Depth density of upper-plate EQ near each cell, ' num2str(startdate) '-' num2str(enddate-1)],[num2str(numperiods) ' temporal samples: wtd. median at ' num2str(depthcens_1km(plotindex)) ' km depth']},'fontsize',20,'fontweight','normal')
     %return
     print(plotindex,['upperplate_' num2str(startdate) 'to' num2str(enddate-1) '_md' num2str(maxdepth) '_timevarMc_' num2str(numperiods) 'periods_minnum' num2str(minnum) '_fix_wmed_' num2str(depthcens_1km(plotindex)) '_finebins_lu' num2str(locuncert) '_' date '.png'],'-dpng','-r280');
end
return
%%
%depthcens = depthcens;
%     plotdepth = depthcens(plotindex);
%     densitygrid_toplot = densitygrid_wmed(:,:,plotindex);
figure(100); clf; hold on; box on;
set(gca,'DataAspectRatio',[1/cosd(-41) 1 111.1/5])
set(gcf,'position',[0 725 525 725])
plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.75 0.75 0.75],'linewidth',0.25)
%plot(eqs_AC.poly50_lons,eqs_AC.poly50_lats,'color','k','linewidth',0.5)
plot(coastlines.X,coastlines.Y,'color',[0 0.25 0.375],'linewidth',0.5)
plot(lakes.X,lakes.Y,'color',[0 0.25 0.375],'linewidth',0.5)

scatter3(hik_CW.lons,hik_CW.lats,-hik_CW.depths,12.5,hik_CW.depths,'v','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
scatter3(puy_HS.lons,puy_HS.lats,-puy_HS.depths,2.5,puy_HS.depths,'^','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
%plot(hik_CW.lons,hik_CW.lats,'*')
%text(hik_CW.lons,hik_CW.lats,num2str([1:length(hik_CW.lats)]'))
pcolor(longrid,latgrid,d90); shading interp;

xlim([165.5 179.5])
ylim([-48 -34])
caxis([0 40])
zlim([-40 0])
colormap(flip(jet4_1000))
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
title({['D90 of upper-plate earthquakes near each cell, ' num2str(startdate) '-' num2str(enddate-1)],['(cross-val. w/ ' num2str(numperiods) ' periods, weighted median)']},'fontsize',20,'fontweight','normal')
%return
print(100,['upperplate_' num2str(startdate) 'to' num2str(enddate-1) '_md' num2str(maxdepth) '_timevarMc_' num2str(numperiods) 'periods_minnum' num2str(minnum) '_fix_wmed_d90.png'],'-dpng','-r280');
%end
return
%%
for plotindex = 2:length(depthcens_1km)
     plotdepth = depthcens_1km(plotindex);
     densitygrid_toplot = depthdensitypercell_w16(:,:,depthcens_1km==plotdepth);
     figure(plotindex); clf; hold on; box on;
     set(gca,'DataAspectRatio',[1/cosd(-41) 1 111.1/5])
     set(gcf,'position',[0 725 530 725])
     %plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.8 0.8 0.8],'linewidth',0.25)
     %plot(eqs_AC.poly50_lons,eqs_AC.poly50_lats,'color','k','linewidth',0.5)
     plot(coastlines(:,1),coastlines(:,2),'color',[0 1/3 1/2],'linewidth',0.5)

     scatter3(hik_CW.lons,hik_CW.lats,-hik_CW.depths,12.5,hik_CW.depths,'v','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
     scatter3(puy_HS.lons,puy_HS.lats,-puy_HS.depths,2.5,puy_HS.depths,'^','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
     %plot(hik_CW.lons,hik_CW.lats,'*')
     %text(hik_CW.lons,hik_CW.lats,num2str([1:length(hik_CW.lats)]'))
     pcolor(longrid,latgrid,densitygrid_toplot); shading interp;

     xlim([165.5 179.5])
     ylim([-48 -34])
     caxis([0 clim_prc])
     zlim([-40 0])
     colormap((jet4w_1000))
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
     ax_width = outerpos(3) - ti(1)*1.33;
     ax_height = outerpos(4) - ti(2) - ti(4);
     ax.Position = [left bottom ax_width ax_height];
     %return
     %return
     title({['Depth density of upper-plate EQ near each cell, ' num2str(startdate) '-' num2str(enddate)],['(cross-val. w/ ' num2str(numperiods) ' periods): ' num2str(plotdepth) ' km depth: weighted 16%']},'fontsize',20,'fontweight','normal')
     %return
     print(plotindex,['./depthdensity_poisson_' num2str(numperiods) 'periods_1yrmin5yr_2ton_histdepths_nonnorm_' num2str(plotdepth) 'km_w16_2_' date '.png'],'-dpng','-r280');
end

%%
for plotindex = 2:length(depthcens_1km)
     plotdepth = depthcens_1km(plotindex);
     densitygrid_toplot = depthdensitypercell_w84(:,:,depthcens_1km==plotdepth);
     figure(plotindex); clf; hold on; box on;
     set(gca,'DataAspectRatio',[1/cosd(-41) 1 111.1/5])
     set(gcf,'position',[0 725 530 725])
     %plot(cfm_1pt0_vertices_NaNs(:,1),cfm_1pt0_vertices_NaNs(:,2),'color',[0.8 0.8 0.8],'linewidth',0.25)
     %plot(eqs_AC.poly50_lons,eqs_AC.poly50_lats,'color','k','linewidth',0.5)
     plot(coastlines(:,1),coastlines(:,2),'color',[0 1/3 1/2],'linewidth',0.5)

     scatter3(hik_CW.lons,hik_CW.lats,-hik_CW.depths,12.5,hik_CW.depths,'v','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
     scatter3(puy_HS.lons,puy_HS.lats,-puy_HS.depths,2.5,puy_HS.depths,'^','linewidth',0.5,'markeredgecolor',[0.6 0.6 0.6]);
     %plot(hik_CW.lons,hik_CW.lats,'*')
     %text(hik_CW.lons,hik_CW.lats,num2str([1:length(hik_CW.lats)]'))
     pcolor(longrid,latgrid,densitygrid_toplot); shading interp;

     xlim([165.5 179.5])
     ylim([-48 -34])
     caxis([0 clim_prc])
     zlim([-40 0])
     colormap((jet4w_1000))
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
     ax_width = outerpos(3) - ti(1)*1.33;
     ax_height = outerpos(4) - ti(2) - ti(4);
     ax.Position = [left bottom ax_width ax_height];
     %return
     %return
     title({['Depth density of upper-plate EQ near each cell, ' num2str(startdate) '-' num2str(enddate)],['(cross-val. w/ ' num2str(numperiods) ' periods): ' num2str(plotdepth) ' km depth: weighted 84%']},'fontsize',20,'fontweight','normal')
     %return
     print(plotindex,['./depthdensity_poisson_' num2str(numperiods) 'periods_1yrmin5yr_2ton_histdepths_nonnorm_' num2str(plotdepth) 'km_w84_2_' date '.png'],'-dpng','-r280');
end

%%
save(['depthdensitygrid_poisson_' num2str(numperiods) 'periods_start2000_1yrmin5yr_2ton_histdepths_nonnorm_' date],'-v7.3')