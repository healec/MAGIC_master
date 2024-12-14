% Profile builder for MAGIC, MAGIC-FOREST, and CGCAM-----------------------
%
% The user can change the options in the "user options to set portion and
% the code will output a profile.data file for MAGIC or MAGIC-FOREST and
% background.wind for CGCAM
close all
clear all
clc

%User options to set-------------------------------------------------------
%--------------------------------------------------------------------------

%Set code flag (1 for MAGIC, 2 for MAGIC-FOREST, 3 for CGCAM)
codeflag=1;

% Choose your altitude range and resolution (km)
altlow=0;%lower altitude in km
althigh=302; %upper altitude in km
dzin=1;  %altitude step size in km
smoothin=3; %smoothing factor in km

%Select the latitude of interest
latin=44.5;
%choose latitude averaging length (+/- from center)
%latext=5;

%Select the longitude of interest
lonin=-100.5;
%choose longitude averaging length (+/- from center)
%lonext=5;

%Select the year of interest
yearin=2022;

%Select the month of interest (1-12)
monthin=7;

%Select the day of interest
dayin=5;

%Select hour of interest (1-24 UT)
hourin=18;

hourarray=hourin;

% Setup MSIS and HWM input arrays------------------------------------------
d=datetime(yearin,monthin,dayin);
d.Format='yyyy-MM-dd'
doyin = day(d,'dayofyear');

H=[altlow:dzin:althigh]'*1000;
YEAR=ones(length(H),1)*yearin;
DOY=ones(length(H),1)*doyin;
LAT=ones(length(H),1)*latin;
LON=ones(length(H),1)*lonin;
SEC=ones(length(H),1)*(3600)*(hourin)/2;
LST=(SEC/3600 + LON/15);

        %Call NRLMSIS and hwm
        [T, RHO] = atmosnrlmsise00(H, LAT, LON, YEAR, DOY, SEC, LST);
        TempMSIS=T(:,2);
        % converts kg/cm^3 to g/cm^3
        dens=RHO(:,6)*(1000/1e6);
        height=H./1000;
        dox=RHO(:,2)./1e6;
        dnit2=RHO(:,3)/1e6;
        dox2=RHO(:,4)/1e6;
        dhyd=RHO(:,7)/1e6;

         heighta=altlow:smoothin:althigh;
         tempa=interp1(height,TempMSIS,heighta,'pchip');
         densa=interp1(height,dens,heighta,'pchip');
         dox2a=interp1(height,dox2,heighta,'pchip');
         doxa=interp1(height,dox,heighta,'pchip');
         dhyda=interp1(height,dhyd,heighta,'pchip');
         dnit2a=interp1(height,dnit2,heighta,'pchip');

          height=(altlow-1.5*dzin):dzin:(althigh+1.5*dzin);

        temp=interp1(heighta,tempa,height,'pchip');
        dens=interp1(heighta,densa,height,'pchip');
        dox2=interp1(heighta,dox2a,height,'pchip');
        dox=interp1(heighta,doxa,height,'pchip');
        dhyd=interp1(heighta,dhyda,height,'pchip');
        dnit2=interp1(heighta,dnit2a,height,'pchip');
counter=0;
for k=1:length(hourarray)

        SEC=ones(length(H),1)*(3600)*hourarray(k);
        W = atmoshwm(LAT, LON, H,'DAY',DOY,'SECONDS',SEC);

        height=H./1000;

        Utemp=squeeze(W(:,2));
        Vtemp=squeeze(W(:,1));

        heighta=altlow:smoothin:althigh;

         winda=interp1(height,Utemp,heighta,'pchip');
         windb=interp1(height,Vtemp,heighta,'pchip');

         height=(altlow-1.5*dzin):dzin:(althigh+1.5*dzin);

         windU=interp1(heighta,winda,height,'pchip');
         windV=interp1(heighta,windb,height,'pchip');

 Na=1e8*exp(-0.5*((height*1e3)-90e3).^2./(3e3^2));
 dne=Na;

 dT_dz=diff(temp)./(dzin*1000);
 N2=(9.81./temp(1:end-1)).*(dT_dz+(9.8/1000));

%make the winds go to zero at the bottom-----------------------------
damper=linspace(0,1,10e3/(dzin*1000));
windU(1:2)=0;
windV(1:2)=0;
windU(3:2+length(damper))=windU(3:2+length(damper)).*damper;
windV(3:2+length(damper))=windV(3:2+length(damper)).*damper;

%final figure
set(figure(2),'Position',[425         958        1940         379])
figure(2)
subplot(1,4,1)
plot(windU,height,windV,height)
xlabel('Wind speed (m/s)')
ylabel('Altitude (km)')
legend('Zonal winds','Meridional winds')
title('Winds (m/s)')
ylim([altlow althigh])

subplot(1,4,2)
plot(temp,height)
xlabel('Temperature (K)')
ylabel('Altitude (km)')
title('Temperature (K)')
ylim([altlow althigh])

subplot(1,4,3)
semilogx(abs(dox),height,'k-',abs(dox2),height,'k-.',abs(dnit2),height,'k:');
legend('O','O_2','N_2');
xlabel('Log Concentrations [Molecules cm^-3]');
ylabel('Altitude [km]');
axis([0 10^20 altlow althigh]);
title('Major Species Concentrations');

subplot(1,4,4)
plot(N2,height(1:end-1))
xlabel('N^2')
ylabel('Altitude (km)')
title('Buoyancy frequency squared (rads/s)^2')
ylim([altlow althigh])

Ufinal(:,counter+1)=windU;
Vfinal(:,counter+1)=windV;

%write file to output
% for MAGIC
if codeflag==1
profile1=[height' dox' dnit2' dox2' dens' temp' dhyd' dne' windU' windV'];
dlmwrite('profile.data',profile1,'delimiter','\t','precision','%.6E','-append');
end
clear windU windV winda windb
counter=counter+1;
end

