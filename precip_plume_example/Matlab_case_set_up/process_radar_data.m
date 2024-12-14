
% Precipitation and Reflectivity Data from MRMS (available from 01/01/2015)-----------------------------------------------------

%insert paramters of interest------------------------------------------


year='2022'; %must be 2015 or later.
month='07'; %Day 186
day='05';
hourstart='16';
minutestart='00'; %must be multiple of 2
hourend='17';
minuteend='40';%must be multiple of 2
dt=1; % choose the array temporal resolution you want (minutes)
dx=2; % choose the array spatial resolution you want (km)
dy=2;
ramptime=40; % choose a ramp up period to avoid transients (mins)
arrayxsize=250; % choose number of grid points for array
arrayysize=200;
loncent=-100.500;%-99.9950; %choose center of array location in lat lon
latcent=44.5;%35;
movie=1; %select 1 if you wish to produce an output movie
moviename='Derecho_5Jul2022_16UT'; %choose a file name for the movie
savefile=1; %select 1 if you wish to save the output as a .mat file
filename='Derecho_5Jul2022_16UT'; %choose a file name for the .mat file

tStart = cputime;

%-------------------------------------------------------------------------

%Create lat lon array-----------------------------------------------------

x=0:dx:dx*arrayxsize;
y=0:dy:0+dy*arrayysize;

% Create new even km latitude array
Lat=latcent-(arrayysize/2)*(1/111)*dy:(1/111)*dy:latcent+(arrayysize/2)*(1/111)*dy;

% Create new even km longitude array
for j=1:length(Lat)
Lon2(:,j)=loncent-(arrayxsize/2)*(1/(cosd(Lat(j)).*111.325))*dx:(1/(cosd(Lat(j)).*111.325))*dx:loncent+(arrayxsize/2)*(1/(cosd(Lat(j)).*111.325))*dx;
end

%rep the new latitude array
Lat2=repmat(Lat,[size(Lon2,1),1]);

tStart = cputime;
%Download data------------------------------------------------------------
%Get path variable for unix system calls
PATH = getenv('PATH');
setenv('PATH', [PATH ':/opt/local/bin/:/usr/local/bin:usr/bin']);

% set path to webserver
s4=strcat('https://mtarchive.geol.iastate.edu/',num2str(year),'/',num2str(month),'/',num2str(day),'/mrms/ncep/PrecipRate/');
s4r=strcat('https://mtarchive.geol.iastate.edu/',num2str(year),'/',num2str(month),'/',num2str(day),'/mrms/ncep/SeamlessHSR/');

%Create a time array between selected dates and times
t1=datetime(str2num(year),str2num(month),str2num(day),str2num(hourstart),str2num(minutestart),0);
t2=datetime(str2num(year),str2num(month),str2num(day),str2num(hourend),str2num(minuteend),0);

%get a list of the files in this directory ---------------
system(['wget --no-check-certificate ',s4])
DD=readlines('index.html');
B=DD(contains(DD,'grib2'));
C=extractBetween(B,"Precip",".gz",'Boundaries','inclusive');
C=C(:,1);
system(['/bin/rm *.html'])

system(['wget --no-check-certificate ',s4r])
DDr=readlines('index.html');
Br=DDr(contains(DDr,'grib2'));
Cr=extractBetween(Br,"Seamless",".gz",'Boundaries','inclusive');
Cr=Cr(:,1);
system(['/bin/rm *.html'])

%Extract just the data and times
Tt=extractBetween(C,"00_",".",'Boundaries','exclusive');
Tt2=datetime(Tt,'InputFormat','yyyyMMdd-HHmmss');

% Find wthe incidies that encompass our time range of interest
[Istart A]=find(t1 <= Tt2);
[Iend A]=find(t2 <= Tt2);
Indstart=Istart(1);
Indend=Iend(1);

%initialize arrays-------------------------------------------
Reflfull=[];
precipfull=[];
timefull=[];

%Perfom time loop to download the files-----------------------------------

for i=Indstart:Indend
%create path to file on webserver
s5=char(Cr(i));
s6=char(C(i));
s=strcat(s4r,s5);
sa=strcat(s4,s6); 

%Get data from webserver
system(['wget --no-check-certificate ',s]);
system(['wget --no-check-certificate ',sa]);
%Unzip the file
system(['/usr/bin/gzip -d ',s5])
system(['/usr/bin/gzip -d ',s6])
%Convert from grib2 to netCDF
system(['wgrib2 ' ,s5(1:end-3),' -netcdf ' ,s5(1:end-8) 'nc']);
system(['wgrib2 ' ,s6(1:end-3),' -netcdf ' ,s6(1:end-8) 'nc']);
%Retrieve variables from netCDF
precip=ncread(strcat(s6(1:end-8),'nc'),'PrecipRate_0mabovemeansealevel');
precip(precip<0)=0;
Refl=ncread(strcat(s5(1:end-8),'nc'),'SeamlessHSR_0mabovemeansealevel');
Refl(Refl<0)=0;
lat=ncread(strcat(s5(1:end-8),'nc'),'latitude');
lon=ncread(strcat(s5(1:end-8),'nc'),'longitude');
time=ncreadatt(strcat(s5(1:end-8),'nc'),'time','reference_date');
system('/bin/rm *.grib2')
system('/bin/rm *.nc')

%Get only the latitude/longitude subset needed-----------------------------

%Find indicies of the array corresponding to the limits
[A Latlow]=min(abs(min(min(Lat2))-lat));
[A Lathigh]=min(abs(max(max(Lat2))-lat));
[A Lonlow]=min(abs(min(min(Lon2))-lon));
[A Lonhigh]=min(abs(max(max(Lon2))-lon));

%define new subsetted lat and lon arrays 
newlat=lat(Latlow:2:Lathigh);
newlon=lon(Lonlow:2:Lonhigh);

%Define subsetted arrays
Reflsub=Refl(Lonlow:2:Lonhigh,Latlow:2:Lathigh);
precipsub=precip(Lonlow:2:Lonhigh,Latlow:2:Lathigh);

% change this
Reflfull=cat(3,Reflfull,Reflsub);
precipfull=cat(3,precipfull,precipsub/6); %divide by 6 to get mm/10 mins
timefull=cat(3,timefull,time);

clear Refl precip Reflsub precipsub

 end


%Interpolate the data to the required size---------------------------------
%Spatial interp from center--------------------------------------------------------

%rep the original lat lon matricies for interpolation
newlatrep=repmat(newlat',[size(newlon,1),1]);
newlonrep=repmat(newlon,[1,size(newlat,1)]);

% Interpolate onto even km grid using griddata.
for k=1:size(Reflfull,3)
precipinterp(:,:,k)=griddata(newlonrep,newlatrep,squeeze(precipfull(:,:,k)),Lon2,Lat2);
Reflinterp(:,:,k)=griddata(newlonrep,newlatrep,squeeze(Reflfull(:,:,k)),Lon2,Lat2);
end

% Interpolate onto even km grid using scatteredInterpolate.
% for k=1:size(Reflfull,3)
%     F = scatteredInterpolant(Lon2,Lat2,squeeze(precipfull(:,:,k)));
%     precipinterp2(:,:,k) = F(newlonrep,newlatrep);
% 
%     FF = scatteredInterpolant(Lon2,Lat2,squeeze(Reflfull(:,:,k)));
%     Reflinterp2(:,:,k) = F(newlonrep,newlatrep);
% end

%Change interpolated NaNs to zeros
Reflinterp(isnan(Reflinterp))=0;
precipinterp(isnan(precipinterp))=0;

%interpolate in time------------------------------------------------------

time_orig=Tt2(Indstart:Indend);
dtorig=minutes(diff(time_orig));
torig=0:(length(time_orig)-1);
for i=2:length(torig)
torig(i)=torig(i-1)+dtorig(i-1);
end

timefinal=0:dt*60:torig(end)*60;
time_real=t1+seconds(timefinal);
%Define ramp up factor
ramp=ones(length(timefinal),1);
ramp(1:ramptime/dt)=linspace(0,1,ramptime/dt);
assert(length(ramp)==length(timefinal),'Error:Ramp time is longer than array')

for i=1:size(precipinterp,1)
    for j=1:size(precipinterp,2)
           precipinterpt(i,j,:) = interp1(torig,squeeze(precipinterp(i,j,:)),timefinal./(dt*60)).*ramp';
           Reflinterpt(i,j,:) = interp1(torig,squeeze(Reflinterp(i,j,:)),timefinal./(dt*60)).*ramp';
    end
end

clear X1 Y1 T1 X2 Y2 T2 Reflinterp precipinterp

 tEnd = cputime - tStart


%Plot the data-------------------------------------------------------------
if movie==1
figure(1)
set(figure(1),'Position',[ 25         209        1391         517])
whitebg('white');
set(figure(1),'DefaultTextFontSize',16)
set(figure(1),'DefaultLineLineWidth',1)
%%writeObj = VideoWriter('simulation.mp4','MPEG-4');
writeObj = VideoWriter(moviename,'MPEG-4');
%
writeObj.FrameRate = 8;
writeObj.Quality = 100;
open(writeObj);

for k=1:size(Reflinterpt,3)
    latlim=[min(min(Lat2)) max(max(Lat2))];
lonlim=[min(min(Lon2)) max(max(Lon2))];


subplot(1,2,1)
ax=worldmap(latlim,lonlim);
mstruct = gcm;
latlim = mstruct.maplatlimit;
lonlim = mstruct.maplonlimit;
geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',1);
geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
 pcolorm(Lat2,Lon2,Reflinterpt(:,:,k))
 shading interp
 shading flat
 caxis([0 70])
 geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',1);
 geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
 colorbar
 xlabel('Longitude')
 ylabel('Latitude')
 title(['Reflectivity (dBz) at ',datestr(time_real(k))])


subplot(1,2,2)
ax=worldmap(latlim,lonlim);
mstruct = gcm;
latlim = mstruct.maplatlimit;
lonlim = mstruct.maplonlimit;
geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',1);
geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
 pcolorm(Lat2,Lon2,precipinterpt(:,:,k))
 shading interp
 shading flat
 caxis([0 25])
 geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',1);
 geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
 colorbar
 xlabel('Longitude')
 ylabel('Latitude')
 title(['Precipitation Rate (mm/10 mins) at ',datestr(time_real(k))])
%xlim([-87.77 -83.74])
%ylim([29 33])
%caxis([0 15])

figure(1)
set(figure(1),'color','w')    
frame = getframe(gcf);
writeVideo(writeObj,frame);
pause(0.09)
end

close(writeObj)
end

%Write the precipitation .txt file for each time step--------------------------------------

for k=1:size(precipinterpt,3)
    strcat('precip',sprintf( '%05d', timefinal(k)),'.txt' )
fileID = fopen(strcat('precip',sprintf( '%05d', timefinal(k)),'.txt' ),'w');
fprintf(fileID,'%2s %2s %2s\n','dxf','dyf','dtf');
fprintf(fileID,'%4u %4u %4u\n',dx*1000,dy*1000,dt*60);
fprintf(fileID,'%4s %4s %2s %2s\n','x','y','precip','precipt+dt');

counter=0;
  for i=1:size(precipinterpt,1)
       for j=1:size(precipinterpt,2)
            if precipinterpt(i,j,k)>0.1
                counter=counter+1;
              if k==size(precipinterpt,3)
                  fprintf(fileID,'%4u %4u %2.1f %2.1f\n',x(i)*1000,y(j)*1000,precipinterpt(i,j,k),0);
              else
                  fprintf(fileID,'%4u %4u %2.1f %2.1f\n',x(i)*1000,y(j)*1000,precipinterpt(i,j,k),precipinterpt(i,j,k+1));
              end
            end
        end
  end
  if counter==0
      fprintf(fileID,'%4u %4u %2.1f %2.1f\n',0,0,0,0);
      fprintf(fileID,'%4u %4u %2.1f %2.1f\n',0,0,0,0);
  end
fclose(fileID);
end



%-Save the .mat file------------------------------------------------------

if savefile==1
save([filename,'.mat'],'precipinterpt','Reflinterpt','Lat2','Lon2','x','y','timefinal','time_real','-v7.3')

end

tEnd = cputime - tStart