
%User inputs---------------------------------------------
%Point to the .mat file produced from the process_radar_data.m script
load('../Derecho_5Jul2022_16UT.mat')
%Point to where the output files are located that you want to plot
load_dir='../Matlab_outputs/Derecho_testcase';
%Name the movie file
moviename='Derecho_testcase_movie';
%Which frame to start
tstart=0;
%Frames to skip
tstep=5;
%Frame to end at
tend=65;

%end of user inputs---------------------------------------------

w = warning ('off','all');
tinputs=tstart:tstep:tend;

    figure(1)
     set(figure(1),'Position',[74         303        1630         963])
     whitebg('white');
     set(figure(1),'DefaultTextFontSize',16)
     set(figure(1),'DefaultLineLineWidth',1)
     writeObj = VideoWriter(moviename,'MPEG-4');
     writeObj.FrameRate = 6;
     writeObj.Quality = 100;
     open(writeObj);

for i=1:length(tinputs)
    Frame=tinputs(i);
    filein=[load_dir,num2str(Frame),'.mat'];
    load(filein)


     latlim=[min(min(Lat)) max(max(Lat))];
     lonlim=[min(min(Lon)) max(max(Lon))];

       subplot(2,2,1)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat2,Lon2,Reflinterpt(:,:,tinputs(i)+1))
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
        title(['Reflectivity (dBz) on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)




      subplot(2,2,2)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,T_out(:,:,1))
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
        title(['T_p (K) at z=',num2str(zref(1)),'km altitude on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)
        if max(max(T_out(:,:,1))) > 0
        clim([-max(max(abs(T_out(:,:,1)))) max(max(abs(T_out(:,:,1))))])
        end

       subplot(2,2,3)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,T_out(:,:,2))
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
        title(['T_p (K) at z=',num2str(zref(2)),'km altitude on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)
        if max(max(T_out(:,:,2))) > 0
        clim([-max(max(abs(T_out(:,:,2)))) max(max(abs(T_out(:,:,2))))])
        end


      subplot(2,2,4)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat,Lon,T_out(:,:,3))
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
        title(['T_p (K) at z=',num2str(zref(3)),'km altitude on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)
         if max(max(T_out(:,:,3))) > 0
        clim([-max(max(abs(T_out(:,:,3)))) max(max(abs(T_out(:,:,3))))])
        end


        figure(1)
        frame = getframe(gcf);
        writeVideo(writeObj,frame);

pause(0.09)

end

close(writeObj)

%Plotting with Xz slices


    figure(2)
     set(figure(2),'Position',[74         825        2456         441])
     whitebg('white');
     set(figure(2),'DefaultTextFontSize',16)
     set(figure(2),'DefaultLineLineWidth',1)
     writeObj = VideoWriter(strcat(moviename,'_xy'),'MPEG-4');
     writeObj.FrameRate = 6;
     writeObj.Quality = 100;
     open(writeObj);

for i=1:length(tinputs)
    Frame=tinputs(i);
    filein=[load_dir,num2str(Frame),'.mat'];
    load(filein)


       subplot(1,3,1)
        ax=worldmap(latlim,lonlim);
        mstruct = gcm;
        latlim = mstruct.maplatlimit;
        lonlim = mstruct.maplonlimit;
        pcolorm(Lat2,Lon2,Reflinterpt(:,:,tinputs(i)+1))
        shading interp
        shading flat
        alpha(0.7)
        caxis([0 65])
        geoshow('usastatehi.shp', 'FaceColor', 'none','LineWidth',0.5);
        geoshow('landareas.shp', 'FaceColor', 'none','LineWidth',0.5);
        colorbar
        hold on
        plotm(ones(size(x))*Lat(1,yref), Lon(:,yref)','-r','LineWidth',2)
        plotm(Lat(xref,:), Lon(xref,yref)*ones(size(y)),'-r','LineWidth',2)
        hold off
        colormap(gca,'jet')
        % Create ylabel
        ylabel('Latitude','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['Reflectivity (dBz) on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)


      subplot(1,3,2)
        imagesc(Lon(:,yref),z./1000,T_outxz')
        axis xy
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Altitude (km)','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Longitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['T_p (K) at Latitude ',num2str(Lat(1,yref)),' on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)
         if max(max(T_outxz)) > 0
        clim([-max(max(abs(T_outxz))) max(max(abs(T_outxz)))])
        end
     

      subplot(1,3,3)
        imagesc(Lat(xref,:),z./1000,T_outyz')
        axis xy
        colormap(gca,gray)
        colorbar
        % Create ylabel
        ylabel('Altitude (km)','LineWidth',1,'FontSize',16);
        % Create xlabel
        xlabel('Latitude','LineWidth',1,'FontSize',16);
        % Create title
        title(['T_p (K) at Longitude ',num2str(Lon(xref,yref)),' on ',datestr(time_real(tinputs(i)+1))],'LineWidth',1,'FontSize',16)
         if max(max(T_outyz)) > 0
        clim([-max(max(abs(T_outyz))) max(max(abs(T_outyz)))])
        end
       


        figure(2)
        frame = getframe(gcf);
        writeVideo(writeObj,frame);

pause(0.09)

end

close(writeObj)


