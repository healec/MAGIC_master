
%----------------------User inputs----------------------------------------
%-------------------------------------------------------------------------

%-----Point to where the fort files are
addpath '../outputs'
%-----Point to the .mat file produced from the process_radar_data.m script
load('../Derecho_5Jul2022_16UT.mat')
%point to where to save the output files and the names to give them
output_dir='../Matlab_outputs/';
output_name='Derecho_testcase';

Frame=0; % Frame number to start
MaxFrames=66; % Number of Frames (same as nout in claw3ez.data)
zref=[35, 55, 87]; %define altitude indexes to extract x-y slices
yref=100; %define y index to extract an x-z slice
xref=50; %define x index to extract an y-z slice

%--------------------End of user inputs------------------------------------
%--------------------------------------------------------------------------

Lat=Lat2(1:end-1,1:end-1);
Lon=Lon2(1:end-1,1:end-1);
x=x(1:end-1);
y=y(1:end-1);

% Import needed slice or make a slice from imported data
filesbasement = 'fort.qq';

for i=1:1:MaxFrames

    nameCur = strcat(filesbasement,num2str(Frame,'%04i'),'id0000','.h5');

    attr = h5readatt(nameCur,'/Pid','Parameters'); % Thread MASTER (0) always outputs its data
    
% New output type
gridno = attr(1);
level = attr(2);
mx = attr(6)/attr(22);
my = attr(7)/attr(23);
mz = attr(8)/attr(24);
dx = attr(18);
dy = attr(19);
dz = attr(20);
t = attr(3);
iframe = attr(25);
meqn = attr(4);
lx = attr(22);
ly = attr(23);
lz = attr(24);
z=0:dz:((lz*mz)-1)*dz;

    fprintf('Working on data for Frame: %d at time: %d s for file: %s\n',Frame,t,nameCur);

    %-------------- Input full 3D domain  --------------%

        cellataltitude = 1;
        id = 1;
        dataset=[];
        datafullset=[];
        for ii=1:1:lx
        for jj=1:1:ly
        nameCur = strcat(filesbasement,num2str(Frame,'%04i'),'id',num2str(id-1,'%04i'),'.h5');
        namedataset = strcat('/Pid',num2str(id));
        tempp = hdf5read(nameCur,'Pid');
       dataset = cat(3,dataset,tempp(:,1:1:end,1:1:end,1:1:end));
    
        id = id+1;
        %id
        end
        datafullset = cat(2,datafullset,dataset);
        size(datafullset)
        dataset=[];
        end
        
datafullset = datafullset(:,1:1:end,1:1:end,:);

fprintf('Full 3D domain is loaded\n');


if (Frame==0)
dox0=squeeze(6.022d23.*datafullset(6,:,:,:).*datafullset(1,:,:,:).*1e-3.*(1/16));
dnit20=squeeze(6.022d23.*(1-datafullset(6,:,:,:)-datafullset(7,:,:,:)).*datafullset(1,:,:,:).*1e-3.*(1/28));
dox20=squeeze(6.022d23.*datafullset(7,:,:,:).*datafullset(1,:,:,:).*1e-3.*(1/32));
gammam=(7.*(dox20+dnit20)+5.*(dox0))./(5.*(dox20+dnit20)+3.*(dox0));
Rm=8.31d3./((1d-6).*((dox0.*16+dnit20.*28+dox20.*32)./(1d-6.*(dox0+dnit20+dox20))));

momnt=datafullset(2:4,:,:,:);
momnt=momnt.*momnt;
kinetic=squeeze(0.5*squeeze(sum(momnt,1))./squeeze(datafullset(1,:,:,:)));
T0=(gammam-1).*((squeeze(datafullset(5,:,:,:))-kinetic)./(squeeze(datafullset(1,:,:,:)).*Rm));

end

rho=squeeze(datafullset(1,:,:,:));
momnt=datafullset(2:4,:,:,:);
momnt=momnt.*momnt;
kinetic=squeeze(0.5*squeeze(sum(momnt,1))./squeeze(datafullset(1,:,:,:)));
Rm=8.31d3./((1d-6).*((dox0.*16+dnit20.*28+dox20.*32)./(1d-6.*(dox0+dnit20+dox20))));
gammam=(7.*(dox20+dnit20)+5.*(dox0))./(5.*(dox20+dnit20)+3.*(dox0));
T=(gammam-1).*((squeeze(datafullset(5,:,:,:))-kinetic)./(squeeze(datafullset(1,:,:,:)).*Rm));
P=rho.*Rm.*T;
       
U=squeeze(datafullset(2,:,:,:)./datafullset(1,:,:,:));
V=squeeze(datafullset(3,:,:,:)./datafullset(1,:,:,:));
W=squeeze(datafullset(4,:,:,:)./datafullset(1,:,:,:));  
dox=(squeeze(6.022d23.*datafullset(6,:,:,:).*datafullset(1,:,:,:).*1e-3.*(1/16)))-dox0;
dnit2=(squeeze(6.022d23.*(1-datafullset(6,:,:,:)-datafullset(7,:,:,:)).*datafullset(1,:,:,:).*1e-3.*(1/28)))-dnit20;
dox2=(squeeze(6.022d23.*datafullset(7,:,:,:).*datafullset(1,:,:,:).*1e-3.*(1/32)))-dox20;
dhyd = squeeze(datafullset(9,:,:,:)).*(dox+dox2+dnit2);
dox3 = squeeze(datafullset(10,:,:,:)).*(dox+dox2+dnit2);
temps=T-T0;  
U_out=U(:,:,zref);
V_out=V(:,:,zref);
W_out=W(:,:,zref);
T_out=temps(:,:,zref);
%P_out=P(:,:,zref);
%rho_out=rho(:,:,zref);
%dnit2_out=dnit2(:,:,zref);
%dox_out=dox(:,:,zref);

U_outxz=squeeze(U(:,yref,:));
V_outxz=squeeze(V(:,yref,:));
W_outxz=squeeze(W(:,yref,:));
T_outxz=squeeze(temps(:,yref,:));

U_outyz=squeeze(U(xref,:,:));
V_outyz=squeeze(V(xref,:,:));
W_outyz=squeeze(W(xref,:,:));
T_outyz=squeeze(temps(xref,:,:));


name = strcat(output_dir,output_name,num2str(Frame),'.mat');
save(name,'U_out','V_out','W_out','T_out','U_outxz','V_outxz','W_outxz','T_outxz','U_outyz','V_outyz','W_outyz','T_outyz','x','y','z','dx','dy','dz','Lat','Lon','zref','yref','xref','-v7.3')

    Frame = Frame + 1;

end