function [PGW,fields_location]=Greedy(fov,n,skymap,resolution,level)
%[PGW,fields_location]=Greedy(fov,n,skymap,resolution,level)
% Greedy finds the locations of the fields using a greedy algorithm;
%
% fov       : the size of field of view, in deg^2. In our model, a square
%             fov is assumed such that the length and the width of the fov
%             is the square root of the fov.
%
% n         : the total number of fields that will be chosen by the greedy
%             algorithm.
%
% skymap    : a string that refers to the name of a txt file. 
%             GW sky localisation information. The skymap should be a
%             healpix file translated into a txt file that contains four columns:
%             1. Index; 2. declination; 3. right ascension; 4. GW prob.
%
% resolution : a number that defines the resolution of the healpix file.
%              The number of pixels of a skymap is equal to 12*2^(2*resolution)
% level : a number that sets the confidence level of the credible region that 
%         will be considered. The default value if 0.9. To use the default
%         value, enter [].
% 
% PGW : a vector containing the values of the sum of the GW probability 
%       within the fields
% 
% fields_location : a n by 2 matrix containing the location infomration of
%                   the selected fields. The first coloumn corresponds to
%                   the declination of the field centers, and the second
%                   the right ascension.

if resolution >10
    [PGW,fields_location]=DownGreedy(fov,n,skymap,resolution,level);
else                                                                        

DelThe = sqrt(fov);                                                         % this sets the boundary within which the pixels of the skymap should be taken into consideration 
DelPhi = sqrt(fov);                                                          

searchrange = sqrt(2) * sqrt(fov)/2 * pi/180;                               % this sets the boundary within which the greedy algorithm can.  


% TotPointNum = 12*2^(2*resolution);

% FOVPointNum = (fov)/(4*180^2/pi)*TotPointNum;                             % the average total points in a field.

message1='Loading the GW information';
disp(message1)
postinfo=load(skymap);                                                      

sorted = sortrows(postinfo,-4);                                             % sorts GW information descending order according to the value of posterior probability.    

maxprobamall=max(sorted(:,4))*(-1);                                                              


textsize=20;

figure
axesm hammer                                                                % defining the axes of the sky map that will be plotted below.
framem; gridm; mlabel; plabel
setm(gca,'MlabelLocation',-150:30:150,'MLineLocation',30,'MlabelParallel',...
    'equator','GLineWidth',.01,'GLineStyle',':','fontsize',textsize,...
    'Plabellocation',-90:15:90)

set(gca,'layer','top')

if isempty(level) == 1
    level=0.9;
elseif isnumeric(level) == 0
    error('level has to either be an empty array of a number.')
end                                                                         % setting the confidence level of the credible region that will be considered

post=sorted(:,4);                                                          

cumulativeSum = cumsum(post) / sum(post);
indexthreshold = find(cumulativeSum >= level, 1, 'first');

message1='% Credible Region of the GW Posterior is being considered, and a sky map is being generated';
display_message=[num2str(level*100),message1];
disp(display_message)
fprintf('\n')

sorted = sorted(1:indexthreshold,:);
dec=90-sorted(:,2)*180/pi;
ra=sorted(:,3)*180/pi;
post=sorted(:,4);

for j=1:length(ra)                                                          % The fits file available to us is written in 0:180 for declination,                                                                                
        if ra(j)>180                                                        % and 0:360 for right ascension, which needs to be in -90:90 for dec, and -180:180 for ra.                                         
            ra(j)=ra(j)-360;  
        end
end
    newcolormap=colormap(flipud(colormap()));
    colormap(newcolormap)
    
    scatterm(dec,ra,40,post,'filled');
    set(gca,'clim',[min(post) max(post)])

    colorbar

[x,y,z] = sph2cart(postinfo(:,3),pi/2-postinfo(:,2),1);                     % Converting the coordinates of all the pixels from spherical coordinate to Cartesian coordinate.

[x1,y1,z1] = sph2cart(sorted(:,3),pi/2-sorted(:,2),1);                      % Converting the coordinates of the pixels within the credible region from spherical coordinate to Cartesian coordinate.


Cartesian = [x,y,z]';                                                      

idx= rangesearch(Cartesian',[x1,y1,z1],searchrange);                        % excluding the pixels that are not within the searchrange of any pixel within the credible region

 num=1;
 PGW=zeros(n,1);

 fields_location=zeros(n,2);

margin=sqrt(4*pi*180^2/pi^2/(12*2^(2*resolution)));
micro_movement1=[0 0 0 -0.5 -0.5 -0.5 0.5 0.5 0.5];
micro_movement2=[-0.5 0 0.5 -0.5 0 0.5 -0.5 0 0.5];

calibration = -0.3;

    if resolution < 10 
         check=zeros(length(idx)*9,1);
        while num<=n                                                                % this while loop looks for the optimal locations for observed fields.
            P_fov=zeros(1,length(idx)*9);                                           % the concept behind it is that every pixel that is within searchrange of any pixel within the credible region will be 
            filtered_idx=cell(1,length(idx)*9);                                     % set as the center of a tentative field. The valus of P_GW for this tentative field is the sum of the GW prob. at all the pixels  
            text1=['Finding the location of the ', num2str(num), 'th field'];       % within the field. after looping over all the relevant pixels, the field rendering the highest value of P_GW will be chosen as the first field.
            disp(text1)                                                             % then the loop will start again with the pixels in the first fields ignored. The process will repeat until n fields are chose.
            filtered_idx_makezero=cell(1,length(idx)*9);
            for i=1:length(idx)

                for z = 1 : 9
                    phi = sorted(i,2)*180/pi+micro_movement1(z)*margin;                                           
                    theta = sorted(i,3)*180/pi+micro_movement2(z)*margin;
                    filtered = postinfo(idx{i},:);
                    filtered_makezeros=postinfo(idx{i},:);

                    newFilter = Cartesian(:,filtered(:,1));
                    newCartUp = roty(phi-DelPhi/2-90)*rotz(theta)*newFilter;
                    newCartLo = roty(phi+DelPhi/2-90)*rotz(theta)*newFilter;

                    newCartLeft = rotz(DelThe/2)*rotx(90-phi)*rotz(theta-90)*newFilter;
                    newCartRight = rotz(-DelThe/2)*rotx(90-phi)*rotz(theta-90)*newFilter;
                    Filter_idx = find(newCartLeft(1,:)>=0 & newCartRight(1,:)<=0 & newCartUp(3,:)<=0 & newCartLo(3,:)>=0);
                    filtered1 = filtered(Filter_idx,:);
                    check(i-1+z)=length(nonzeros(filtered1(:,4)));
                    % filter to get points with appropriate phi

                    P_fov(i-1+z) = sum(filtered1(:,4));

                    newCartUp_makezeros = roty(phi-DelPhi/2-calibration*margin-90)*rotz(theta)*newFilter;
                    newCartLo_makezeros = roty(phi+DelPhi/2+calibration*margin-90)*rotz(theta)*newFilter;

                    newCartLeft_makezeros = rotz(DelThe/2+calibration*margin)*rotx(90-phi)*rotz(theta-90)*newFilter;
                    newCartRight_makezeros= rotz(-DelThe/2-calibration*margin)*rotx(90-phi)*rotz(theta-90)*newFilter;

                    Filter_idx_makezeros = find(newCartLeft_makezeros(1,:)>=0 & newCartRight_makezeros(1,:)<=0 & newCartUp_makezeros(3,:)<=0 & newCartLo_makezeros(3,:)>=0); 

                    filtered_makezeros=filtered_makezeros(Filter_idx_makezeros,:);
                    filtered_idx{i-1+z}=filtered1(:,1);
                    filtered_idx_makezero{i-1+z}=filtered_makezeros(:,1);
                end
            end

            [~,max_id]=max(P_fov);

            postinfo(filtered_idx_makezero{max_id},4)=0;


            declowlim=90-max(postinfo(filtered_idx{max_id},2))*180/pi;
            decuplim=90-min(postinfo(filtered_idx{max_id},2))*180/pi;


            ra_within=postinfo(filtered_idx{max_id},3)*180/pi;

            for oo=1:length(ra_within)
                if ra_within(oo)>180
                    ra_within(oo)=ra_within(oo)-360;
                end
            end



            ra_uplim=max(ra_within);
            ra_lowlim=min(ra_within);

            PGW(num)=P_fov(max_id);

            text1=['The number of non-zero pixels in this field is:',...            % This is checking the number of pixels in a field is of the right other, because the number may not always be the same due to fluctuation.
                num2str(check(max_id))]; 
            disp(text1)


            dec_location=sorted(max_id,2)*180/pi;
            ra_location=sorted(max_id,3)*180/pi;
            dec_location=90-dec_location;
            if ra_location>180
                ra_location=ra_location-360;
            end

            fields_location(num,1)=dec_location;
            fields_location(num,2)=ra_location;

            message2='The location of this field is at ';
            message3=['(dec = ',num2str(dec_location),',',' ra = ', num2str(ra_location),')'];
            disp(message2)
            disp(message3)
            fprintf('\n')


            if ra_uplim-ra_lowlim<300
                h=surfacem([declowlim,decuplim],[ra_lowlim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
                h=surfacem([decuplim,decuplim],[ra_lowlim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
                h=surfacem([declowlim,decuplim],[ra_uplim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
                h=surfacem([declowlim,declowlim],[ra_lowlim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            else
                h=surfacem([declowlim,decuplim],[ra_lowlim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
                h=surfacem([decuplim,decuplim],[ra_uplim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
                h=surfacem([declowlim,decuplim],[ra_uplim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
                h=surfacem([declowlim,declowlim],[ra_uplim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            end
            textm(dec_location,ra_location,num2str(num),'color','k','fontsize',10)
            num=num+1;
        end

    else
         check=zeros(length(idx),1);
        while num<=n                                                                % this while loop looks for the optimal locations for observed fields.
        P_fov=zeros(1,length(idx));                                             % the concept behind it is that every pixel that is within searchrange of any pixel within the credible region will be 
        filtered_idx=cell(1,length(idx));                                       % set as the center of a tentative field. The valus of P_GW for this tentative field is the sum of the GW prob. at all the pixels  
        text1=['Finding the location of the ', num2str(num), 'th field'];       % within the field. after looping over all the relevant pixels, the field rendering the highest value of P_GW will be chosen as the first field.
        disp(text1)                                                             % then the loop will start again with the pixels in the first fields ignored. The process will repeat until n fields are chose.
        for i=1:length(idx)
            phi = sorted(i,2)*180/pi;                                           
            theta = sorted(i,3)*180/pi;
            filtered = postinfo(idx{i},:);
            filtered_makezeros=postinfo(idx{i},:);

            newFilter = Cartesian(:,filtered(:,1));
            newCartUp = roty(phi-DelPhi/2-90)*rotz(theta)*newFilter;
            newCartLo = roty(phi+DelPhi/2-90)*rotz(theta)*newFilter;

            newCartLeft = rotz(DelThe/2)*rotx(90-phi)*rotz(theta-90)*newFilter;
            newCartRight = rotz(-DelThe/2)*rotx(90-phi)*rotz(theta-90)*newFilter;
            Filter_idx = find(newCartLeft(1,:)>=0 & newCartRight(1,:)<=0 & newCartUp(3,:)<=0 & newCartLo(3,:)>=0);
            filtered1 = filtered(Filter_idx,:);
            check(i)=length(nonzeros(filtered1(:,4)));
            % filter to get points with appropriate phi

            P_fov(i) = sum(filtered1(:,4));

            newCartUp_makezeros = roty(phi-DelPhi/2-calibration*margin-90)*rotz(theta)*newFilter;
            newCartLo_makezeros = roty(phi+DelPhi/2+calibration*margin-90)*rotz(theta)*newFilter;

            newCartLeft_makezeros = rotz(DelThe/2+calibration*margin)*rotx(90-phi)*rotz(theta-90)*newFilter;
            newCartRight_makezeros= rotz(-DelThe/2-calibration*margin)*rotx(90-phi)*rotz(theta-90)*newFilter;

            Filter_idx_makezeros = find(newCartLeft_makezeros(1,:)>=0 & newCartRight_makezeros(1,:)<=0 & newCartUp_makezeros(3,:)<=0 & newCartLo_makezeros(3,:)>=0); 

            filtered_makezeros=filtered_makezeros(Filter_idx_makezeros,:);        
            filtered_idx{i}=filtered1(:,1);
            filtered_idx_makezero{i}=filtered_makezeros(:,1);
        end

        [~,max_id]=max(P_fov);
        postinfo(filtered_idx_makezero{max_id},4)=0;


        declowlim=90-max(postinfo(filtered_idx{max_id},2))*180/pi;
        decuplim=90-min(postinfo(filtered_idx{max_id},2))*180/pi;


        ra_within=postinfo(filtered_idx{max_id},3)*180/pi;

        for oo=1:length(ra_within)
            if ra_within(oo)>180
                ra_within(oo)=ra_within(oo)-360;
            end
        end

        ra_uplim=max(ra_within);
        ra_lowlim=min(ra_within);

        PGW(num)=P_fov(max_id);

        text1=['The number of non-zero pixels in this field is:',...            % This is checking the number of pixels in a field is of the right other, because the number may not always be the same due to fluctuation.
            num2str(check(max_id))]; 
        disp(text1)

        dec_location=sorted(max_id,2)*180/pi;
        ra_location=sorted(max_id,3)*180/pi;
        dec_location=90-dec_location;
        if ra_location>180
            ra_location=ra_location-360;
        end


        fields_location(num,1)=dec_location;
        fields_location(num,2)=ra_location;

        message2='The location of this field is at ';
        message3=['(dec = ',num2str(dec_location),',',' ra = ', num2str(ra_location),')'];
        disp(message2)
        disp(message3)
        fprintf('\n')


        if ra_uplim-ra_lowlim<300
            h=surfacem([declowlim,decuplim],[ra_lowlim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            h=surfacem([decuplim,decuplim],[ra_lowlim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            h=surfacem([declowlim,decuplim],[ra_uplim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            h=surfacem([declowlim,declowlim],[ra_lowlim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
        else
            h=surfacem([declowlim,decuplim],[ra_lowlim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            h=surfacem([decuplim,decuplim],[ra_uplim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            h=surfacem([declowlim,decuplim],[ra_uplim,ra_uplim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
            h=surfacem([declowlim,declowlim],[ra_uplim,ra_lowlim],maxprobamall,'linestyle','-','linewidth',2,'facecolor','k');
        end
        textm(dec_location,ra_location,num2str(num),'color','k','fontsize',10)
        num=num+1;
        end
    end
end
