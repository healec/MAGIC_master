str1='Derecho_Output_attempt2_timedep';
step=5
tend=300


load('/Users/healec/Dropbox (Personal)/Convective_sims/Derecho_5th_July_2022/RADAR_from_16UT/Derecho_5Jul2022_16UT.mat')

for i=1:step:tend
    Frame=(i-1)*60;
    filein=[str1,num2str(Frame),'.mat'];
    load(filein,'U_outxz')
    figure(1)
    plot(0:499,U_outxz(1,:))
    hold on

end


filein=[str1,num2str(0),'.mat'];
load(filein,'T_out')
T_out0=T_out;


    figure(1)
     set(figure(1),'Position',[74         303        1630         963])
     whitebg('white');
     set(figure(1),'DefaultTextFontSize',16)
     set(figure(1),'DefaultLineLineWidth',1)
     writeObj = VideoWriter('Derecho_movie_attempt2timedep','MPEG-4');
     writeObj.FrameRate = 6;
     writeObj.Quality = 100;
     open(writeObj);

for i=1:step:tend
    Frame=(i-1)*60;
    filein=[str1,num2str(Frame),'.mat'];
    load(filein)



     latlim=[min(min(Lat)) max(max(Lat))];
     lonlim=[min(min(Lon)) max(max(Lon))];

       subplot(2,2,1)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat2,Lon2,Reflinterpt(:,:,i))
        shading interp
        shading flat
        alpha(0.7)
        caxis([0 65])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colorbar
        colormap(gca,'jet')
        % Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['Reflectivity (dBz) on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)




      subplot(2,2,2)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,T_out(:,:,3)-T_out0(:,:,3))
        shading interp
        shading flat
      %  caxis([155-35 155+35])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['Temperature (K) at z=87km altitude on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)
        clim([-40 40])

       subplot(2,2,3)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,T_out(:,:,5)-T_out0(:,:,5))
        shading interp
        shading flat
      %  caxis([520-200 520+200])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['Temperature (K) at z=130km altitude on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)
        clim([-170 170])


      subplot(2,2,4)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,T_out(:,:,end)-T_out0(:,:,end))
        shading interp
        shading flat
     %   caxis([1090-150 1090+150])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['Temperature (K) at z=250km altitude on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)
        clim([-100 100])


        figure(1)
        frame = getframe(gcf);
        writeVideo(writeObj,frame);

pause(0.09)

end

close(writeObj)

%Plotting with Xz slices

%str1='Derecho_Output';
%str2='Derecho_OutputXZ45';
%load('/Users/healec/Dropbox (Personal)/Convective_sims/Derecho_5th_July_2022/RADAR/Derecho_5Jul2022.mat')
%Lat=Lat2;
%Lon=Lon2;

    figure(2)
     set(figure(2),'Position',[74         503        1734         763])
     whitebg('white');
     set(figure(2),'DefaultTextFontSize',16)
     set(figure(2),'DefaultLineLineWidth',1)
     writeObj = VideoWriter('Derecho_movie_WXZ45_zoom_timedep','MPEG-4');
     writeObj.FrameRate = 8;
     writeObj.Quality = 100;
     open(writeObj);

for i=1:step:tend
    Frame=(i-1)*60;
    filein=[str1,num2str(Frame),'.mat'];
    load(filein,'W_out')

        filein2=[str1,num2str(Frame),'.mat'];
    load(filein2,'W_outxz')

latin=450;

     latlim=[42.5 48];
     lonlim=[-106 -90];

       subplot(2,2,1)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat2,Lon2,Reflinterpt(:,:,i))
        shading interp
        shading flat
        alpha(0.7)
        caxis([0 65])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colorbar
        hold on
        plotm(ones(size(x))*Lat(1,latin), Lon(:,latin)','-r','LineWidth',2)
        hold off
        colormap(gca,'jet')
        % Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['Reflectivity (dBz) on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)


       subplot(2,2,3)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,W_out(:,:,5))
        shading interp
        shading flat
        caxis([-150 150])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colormap(gca,gray)
        colorbar
        hold on
        plotm(ones(size(x))*Lat(1,latin), Lon(:,latin)','-r','LineWidth',2)
        hold off
%         Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
%         Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
%         Create title
        title(['W (m/s) at z=130km altitude on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)


      subplot(2,2,2)
        imagesc(Lon(:,latin),120:350,W_outxz(:,120:350)')
        axis xy
        clim([-300 300])
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Altitude (km)','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['W (m/s) at Latitude ',num2str(Lat(1,latin)),' on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)
        xlim([-106 -90])

      subplot(2,2,4)
        imagesc(Lon(:,latin),0:125,W_outxz(:,1:125)')
        axis xy
        clim([-80 80])
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Altitude (km)','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['W (m/s) at Latitude ',num2str(Lat(1,latin)),' on ',datestr(time_real(i))],'LineWidth',1,'FontSize',16)
        xlim([-106 -90])


        figure(2)
        frame = getframe(gcf);
        writeVideo(writeObj,frame);

pause(0.09)

end

close(writeObj)


%Wavelet analysis---------------------------------------------------------

%for diagonal line---------------------------------------------------------

% figure(8)
% set(figure(8),'Position',[81         960        2456         367])
%  whitebg('white');
%      set(figure(8),'DefaultTextFontSize',16)
%      set(figure(8),'DefaultLineLineWidth',1)
%      writeObj = VideoWriter('Apr_13_2020_movieT250km_Wavelet','MPEG-4');
%      writeObj.FrameRate = 8;
%      writeObj.Quality = 100;
%      open(writeObj);
%
%
% for i=90:408
%     Frame=(i-1)*60;
%     filein=[str1,num2str(Frame),'.mat'];
%     load(filein)
%
%
% for i=1:1200
% Tp_line(i)=squeeze(T(i,i,1));
% Lat_Line(i)=squeeze(Lat(i,i));
% Lon_Line(i)=squeeze(Lon(i,i));
% end
%
%
% dxline=sqrt(dx.^2+dx.^2);
% xline=0:dxline:(length(Tp_line)-1)*dxline;
%
% [wave4,Period4,scale4,COI4]=wavelet2(Tp_line,squeeze(dxline),1,0.01,-1,-1,-1,-1);
% New_Period4=linspace(min(Period4), max(Period4), length(Period4));
% Period4a=repmat(Period4',[1,size(wave4,2)]);
%
% [t2d2, Period2d2]=meshgrid(xline,Period4);
% [New_t2d2, New_Period2d2]=meshgrid(xline,New_Period4);
% waveinterp4=interp2(t2d2, Period2d2,wave4,New_t2d2, New_Period2d2);
%
%
% wave_tot(:,:,i)=abs(waveinterp4).^2./Period4a;
%
%
% subplot(1,3,1)
%      latlim=[min(min(Lat)) max(max(Lat))];
%      lonlim=[min(min(Lon)) max(max(Lon))];
%         ax=worldmap(latlim,lonlim);
%         mstruct = gcm;
%         latlim = mstruct.maplatlimit;
%         lonlim = mstruct.maplonlimit;
%         pcolorm(Lat,Lon,T(:,:,1))
%         shading interp
%         shading flat
%         caxis([125 225])
%         caxis([930-300 930+300])
%         geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
%         geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
%         colormap(gray)
%         colorbar
%         % Create ylabel
%         ylabel('Latitude','LineWidth',1,'FontSize',16);
%         % Create xlabel
%         xlabel('Longitude','LineWidth',1,'FontSize',16);
%         % Create title
%         title(['Temperature (K) at z=90km altitude and ',datestr(Frame/24/60/60,'HH:MM:SS'),' UT'],'LineWidth',1,'FontSize',16)
%         hold on
% plotm(Lat_Line,Lon_Line,'-r','Linewidth',1)
% hold off
% colorbar
% xlabel('Longitude')
% ylabel('Latitude')
% title(['Temperature (K) at z=90km altitude and ',datestr(Frame/24/60/60,'HH:MM:SS'),' UT'],'LineWidth',1,'FontSize',16)
% title(['Temperature (K) at z=250 km altitude and ',datestr(Frame/24/60/60,'HH:MM:SS'),' UT'],'LineWidth',1,'FontSize',16)
%
% colormap(gca,gray)
%
%
% subplot(1,3,2)
% plot(Lon_Line,Tp_line,'Linewidth',1)
% xlabel('Longitude along red line')
% ylabel('T (K)')
% title('Amplitude along red line')
% xlim([min(Lon_Line) max(Lon_Line)])
% colorbar
% ylim([125 225])
% ylim([930-300 930+300])
%
% subplot(1,3,3)
% imagesc(Lon_Line,New_Period4./1000,abs(waveinterp4).^2./Period4a);
% axis xy
% xlabel('Longitude along red line')
% ylabel('Horizontal wavelength (km)')
% title(['Wavelet power spectra along red line'])
% ylim([0 1000])
% %caxis([0 0.3])
% xlim([min(Lon_Line) max(Lon_Line)])
% colorbar
%
%
%        figure(8)
%         frame = getframe(gcf);
%         writeVideo(writeObj,frame);
%
% pause(0.09)
% end
%
% close(writeObj)
%
% figure(9)
% imagesc(Lon_Line,New_Period4./1000,mean(wave_tot,3));
% axis xy
% xlabel('Longitude along red line')
% ylabel('Horizontal wavelength (km)')
% title(['Time averaged wavelet power spectra along red line'])
% ylim([0 1000])
% xlim([min(Lon_Line) max(Lon_Line)])
% colorbar
