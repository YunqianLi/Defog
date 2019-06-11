% Date:20/6/2013
% This code was writen by Changkai Zhao
% Email:changkaizhao1006@gmail.com
% completed and corrected in 20/06/2013 and in 23/06/2013
%
% This algorithm is described in details in 
%
% "Single Image Haze Removal Using Dark Channel Prior",
% by Kaiming He Jian Sun Xiaoou Tang,
% In: CVPR 2009

% details about guilded image filter in
% "Guilded image filtering"
% by Kaiming He Jian Sun Xiaoou Tang

% OUTPUT:
% J is the obtained clear image after visibility restoration
% tmap is the raw transmission map 
% tmap_ref is the refined transmission map .
function [ J,tmap,tmap_ref ] = darkChannel( I,px,w)
close all 
tic;
if (nargin < 3)
   w = 0.95;        % by default, constant parameter w is set to 0.95.
end
if (nargin < 2)
   px = 15;         % by default, the window size is set to 15.
end
if (nargin < 1)
	msg1 = sprintf('%s: Not input.', upper(mfilename));
        eid = sprintf('%s:NoInputArgument',mfilename);
        error(eid,'%s %s',msg1);
end
if ((w >= 1.0) || (w<=0.0))
	msg1 = sprintf('%s: w is out of bound.', upper(mfilename));
        msg2 = 'It must be an float between 0.0 and 1.0';
        eid = sprintf('%s:outOfRangeP',mfilename);
        error(eid,'%s %s',msg1,msg2);
end

if (px < 1)
	msg1 = sprintf('%s: px is out of bound.', upper(mfilename));
        msg2 = 'It must be an integer higher or equal to 1.';
        eid = sprintf('%s:outOfRangeSV',mfilename);
        error(eid,'%s %s',msg1,msg2);
end


% Pick the top 0.1% brightest pixels in the dark channel.
Im=im2double(I);
[dimr,dimc,col]=size(I);
dx=floor(px/2);
% Initial three matrices
J=zeros(dimr,dimc,col);
t_map=zeros(dimr,dimc);
J_darktemp=zeros(dimr,dimc);
tmap_ref=zeros(dimr,dimc);
    if(col==3)
        A_r=0;
        A_g=0;
        A_b=0;
        %% Estimate the atmospheric light (color)
        J_darkchannel=min(Im,[],3);
        for i=(1:dimr)
            for j=(1:dimc)
                ilow=i-dx;ihigh=i+dx;
                jlow=j-dx;jhigh=j+dx;
                if(i-dx<1)
                    ilow=1;
                end
                if(i+dx>dimr)
                    ihigh=dimr;
                end
                if(j-dx<1)
                    jlow=1;
                end
                if(j+dx>dimc)
                    jhigh=dimc;
                end
                J_darktemp(i,j)= min(min(J_darkchannel(ilow:ihigh,jlow:jhigh)));
            end
        end
        J_darkchannel=J_darktemp;
               %get the 0.1% most brightest pixels in the darkchannel
                lmt = quantile(J_darkchannel(:),[.999]);      
                [rind,cind]=find(J_darkchannel>=lmt);
                [enum,~]=size(rind);
                for i=(1:enum)
                    A_r=Im(rind(i),cind(i),1)+A_r;
                    A_g=Im(rind(i),cind(i),2)+A_g;
                    A_b=Im(rind(i),cind(i),3)+A_b;
                end
                %here we get the airlight
                Airlight=[A_r/enum,A_g/enum,A_b/enum];

        %% Estimating the raw transmission map(color)
        Im_n(:,:,1)=Im(:,:,1)./Airlight(1);
        Im_n(:,:,2)=Im(:,:,2)./Airlight(2);
        Im_n(:,:,3)=Im(:,:,3)./Airlight(3);
        tmap=min(Im_n,[],3);
        for i=(1:dimr)
            for j=(1:dimc)
                ilow=i-dx;ihigh=i+dx;
                jlow=j-dx;jhigh=j+dx;
                if(i-dx<1)
                    ilow=1;
                end
                if(i+dx>dimr)
                    ihigh=dimr;
                end
                if(j-dx<1)
                    jlow=1;
                end
                if(j+dx>dimc)
                    jhigh=dimc;
                end
                t_map(i,j)= 1-w*min(min(tmap(ilow:ihigh,jlow:jhigh)));
            end
        end
        
        %% Refine the raw transmission map(color)
        % The auther used soft mapping method to refine t map originally, afterward switch to 
        % guilded filter to refine t map, here we use guilded filter to refine t map.
       
        %tmap_ref=guidedfilter_color(Im,t_map,40,0.001);

        % Alternatively,we can use softmatting to refine the raw t map
        tmap_ref=softmatting(Im,t_map);
        %% Getting the clear image(color)
        [rind,cind]=find(tmap_ref<0.1);
        [enum,~]=size(rind);
        for i=(1:enum)
            tmap_ref(rind(i),cind(i))=0.1;
        end
        J(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_ref+Airlight(1);
        J(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_ref+Airlight(2);
        J(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_ref+Airlight(3);
        
        tmap_ref1 = max(t_map, 0.1);
        J1(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_ref1+Airlight(1);
        J1(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_ref1+Airlight(2);
        J1(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_ref1+Airlight(3);

        figure,imshow(Im),title('Input Image');
        figure,imshow(t_map),title('Raw t map');imwrite(t_map,'t_map.tif');       
        figure,imshow(tmap_ref),title('Refined t map');imwrite(tmap_ref,'tmap_ref.tif');
        figure,imshow(tmap_ref1),title('Not Refined t map');imwrite(tmap_ref1,'tmap_ref1.tif');        
        figure,imshow(J),title('Output Image with soft matting');
        figure,imshow(J1),title('Output Image without soft matting');
        

    end
            psnr(J,J1)
    if(col==1)
        Airlight=0;
        %% Estimate the atmospheric light (gray)
        for i=(1:dimr)
            for j=(1:dimc)
                ilow=i-dx;ihigh=i+dx;
                jlow=j-dx;jhigh=j+dx;
                if(i-dx<1)
                    ilow=1;
                end
                if(i+dx>dimr)
                    ihigh=dimr;
                end
                if(j-dx<1)
                    jlow=1;
                end
                if(j+dx>dimc)
                    jhigh=dimc;
                end
                J_darkchannel(i,j)= min(min(Im(ilow:ihigh,jlow:jhigh)));                                      
            end
        end
               %get the 0.1% most brightest pixels in the darkchannel
                lmt = quantile(J_darkchannel(:),[.999]);      
                [rind,cind]=find(J_darkchannel>=lmt);
                [enum,~]=size(rind);
                for i=(1:enum)
                    Airlight=Im(rind(i),cind(i))+Airlight;
                end
                %here we get the airlight
                Airlight=Airlight/enum;

        %% Estimating the raw transmission map(gray)
        %Im_n=Im./Airlight;

        for i=(1:dimr)
            for j=(1:dimc)
%                 ilow=i-dx;ihigh=i+dx;
%                 jlow=j-dx;jhigh=j+dx;
%                 if(i-dx<1)
%                     ilow=1;
%                 end
%                 if(i+dx>dimr)
%                     ihigh=dimr;
%                 end
%                 if(j-dx<1)
%                     jlow=1;
%                 end
%                 if(j+dx>dimc)
%                     jhigh=dimc;
%                 end
%                 tmap(i,j)= 1-w*min(min(Im_n(ilow:ihigh,jlow:jhigh)));
                  tmap(i,j)= 1-w*(J_darkchannel(i,j)/Airlight);
            end
        end

        
        %% Refine the raw transmission map(gray)
        % The auther used soft mapping method to refine t map originally, afterward switch to 
        % guilded filter to refine t map, here we use guilded filter to refine t map.
        tmap_ref=guidedfilter(Im,tmap,20,0.01);
        %% Getting the clear image(gray)
        [rind,cind]=find(tmap_ref<0.1);
        [enum,~]=size(rind);
            for i=(1:enum)
                tmap_ref(rind(i),cind(i))=0.1;
            end
        J= (Im-Airlight)./tmap_ref+Airlight;


        figure,imshow(Im),title('Input Image');
        figure,imshow(tmap),title('Raw t map');
        figure,imshow(tmap_ref),title('Refined t map');
        figure,imshow(J),title('Output Image');
    end
   toc;
end



