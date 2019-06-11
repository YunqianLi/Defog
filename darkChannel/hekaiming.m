clear
close all

tic;

%% init
w = 0.75;            % by default, constant parameter w is set to 0.95.
windowsize = 15;     % by default, the window size is set to 15.
t0 = 0.1;            % by default, the minimun transmission value
position = [345,180,30,30];

%% load image
I = imread('./images/nanjing.tif');
I = imresize(I, 0.5) * 1.5;
Im=im2double(I);
[rows,cols,channels]=size(I);
figure, imshow(Im),title('origianl image'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(Im,'Im.tif');

%% compute min/darkchannel in pixels
I_darkchannel_pixels=min(Im,[],3);
figure, imshow(I_darkchannel_pixels),title('darkchannel in pixels');
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(I_darkchannel_pixels,'I_darkchannel_pixels.tif');

%% compute min/darkchannel in windows
I_darkchannel_temp=zeros(rows,cols);
dx=floor(windowsize/2);
for i=(1:rows)
    for j=(1:cols)
        ilow=i-dx;ihigh=i+dx;
        jlow=j-dx;jhigh=j+dx;
        if(i-dx<1)
            ilow=1;
        end
        if(i+dx>rows)
            ihigh=rows;
        end
        if(j-dx<1)
            jlow=1;
        end
        if(j+dx>cols)
            jhigh=cols;
        end
        I_darkchannel_temp(i,j)= min(min(I_darkchannel_pixels(ilow:ihigh,jlow:jhigh)));
    end
end
I_darkchannel_windows=I_darkchannel_temp;
figure, imshow(I_darkchannel_windows),title('darkchannel in windows'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(I_darkchannel_windows,'I_darkchannel_windows.tif');


%% Estimate the atmospheric light (color)
% Find the top 0.1% brightest pixels' locations in the dark channel.
A_r=0;
A_g=0;
A_b=0;
lmt = quantile(I_darkchannel_windows(:),[.999]);
[rind,cind]=find(I_darkchannel_windows>=lmt);
[enum,~]=size(rind);
for i=(1:enum)
    A_r=Im(rind(i),cind(i),1)+A_r;
    A_g=Im(rind(i),cind(i),2)+A_g;
    A_b=Im(rind(i),cind(i),3)+A_b;
end
Airlight=[A_r/enum,A_g/enum,A_b/enum];
% Airlight = [0.75,0.75,0.75];
%% Estimating the raw transmission map(color) in pixels
Im_n(:,:,1)=Im(:,:,1)./Airlight(1);
Im_n(:,:,2)=Im(:,:,2)./Airlight(2);
Im_n(:,:,3)=Im(:,:,3)./Airlight(3);
tmap_pixels_temp=min(Im_n,[],3);
tmap_pixels = 1-w*tmap_pixels_temp;
tmap_pixels = max(tmap_pixels, t0);
figure, imshow(tmap_pixels),title('transmission map in pixels'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(tmap_pixels,'tmap_pixels.tif');

%% Bilateral filter the raw transmission map in pixels
r = 5;% 滤波半径
a = 1;% 全局方差
b = 1;% 局部方差
tmap_pixels_filtered = my_bfilt_gray(tmap_pixels,r,a,b);
figure, imshow(tmap_pixels_filtered),title('filtered transmission map in pixels'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(tmap_pixels_filtered,'tmap_pixels_filter.tif');

%% Estimating the raw transmission map(color) in windows
tmap_windows = zeros(rows,cols);
for i=(1:rows)
    for j=(1:cols)
        ilow=i-dx;ihigh=i+dx;
        jlow=j-dx;jhigh=j+dx;
        if(i-dx<1)
            ilow=1;
        end
        if(i+dx>rows)
            ihigh=rows;
        end
        if(j-dx<1)
            jlow=1;
        end
        if(j+dx>cols)
            jhigh=cols;
        end
        tmap_windows(i,j)= 1-w*min(min(tmap_pixels_temp(ilow:ihigh,jlow:jhigh)));
    end
end
tmap_windows = max(tmap_windows,t0);
figure, imshow(tmap_windows),title('transmission map in windows'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(tmap_windows,'tmap_windows.tif');

%% Bilateral filter the raw transmission map in windows
tmap_windows_filtered = my_bfilt_gray(tmap_windows,r,a,b);
figure, imshow(tmap_windows_filtered),title('filtered transmission map in windows'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(tmap_windows_filtered,'tmap_windows_filter.tif');

%% Refine the raw transmission map(color)
% The auther used soft mapping method to refine t map originally, afterward switch to
% guilded filter to refine t map, here we use guilded filter to refine t map.
%tmap_ref=guidedfilter_color(Im,t_map,40,0.001);
% Alternatively,we can use softmatting to refine the raw t map
tmap_refined=softmatting(Im,tmap_windows);
tmap_refined = max(tmap_refined,t0);
figure, imshow(tmap_refined),title('refined transmission map in windows'); 
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(tmap_refined,'tmap_ref.tif');

%% Recover the clear image(color)
% use tmap_pixels
J_pixels=zeros(rows,cols,channels);
J_pixels(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_pixels+Airlight(1);
J_pixels(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_pixels+Airlight(2);
J_pixels(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_pixels+Airlight(3);
figure,imshow(J_pixels),title('recoverd image with tmap\_pixels');
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
imwrite(J_pixels, 'J_pixels.tif');
% use tmap_pixels_filtered
J_pixels_filtered=zeros(rows,cols,channels);
J_pixels_filtered(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_pixels_filtered+Airlight(1);
J_pixels_filtered(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_pixels_filtered+Airlight(2);
J_pixels_filtered(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_pixels_filtered+Airlight(3);
figure,imshow(J_pixels_filtered),title('recoverd image with tmap\_pixels\_filtered');
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
PSNR1 = psnr(J_pixels,J_pixels_filtered);
imwrite(J_pixels_filtered, 'J_pixels_filtered.tif');
% use tmap_windows
J_windows=zeros(rows,cols,channels);
J_windows(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_windows+Airlight(1);
J_windows(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_windows+Airlight(2);
J_windows(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_windows+Airlight(3);
figure,imshow(J_windows),title('recoverd image with tmap\_windows');
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
PSNR2 = psnr(J_pixels,J_windows);
imwrite(J_windows, 'J_windows.tif');
% use tmap_windows_filtered
J_windows_filtered=zeros(rows,cols,channels);
J_windows_filtered(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_windows_filtered+Airlight(1);
J_windows_filtered(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_windows_filtered+Airlight(2);
J_windows_filtered(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_windows_filtered+Airlight(3);
figure,imshow(J_windows_filtered),title('recoverd image with tmap\_windows\_filtered');
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
PSNR3 = psnr(J_windows,J_windows_filtered);
imwrite(J_windows_filtered, 'J_windows_filtered.tif');
% use tmap_refined
J_refined=zeros(rows,cols,channels);
J_refined(:,:,1) = (Im(:,:,1)-Airlight(1))./tmap_refined+Airlight(1);
J_refined(:,:,2) = (Im(:,:,2)-Airlight(2))./tmap_refined+Airlight(2);
J_refined(:,:,3) = (Im(:,:,3)-Airlight(3))./tmap_refined+Airlight(3);
figure,imshow(J_refined),title('recoverd image with tmap\_refined');
hold on;
rectangle('Position',position,'linewidth',3,'EdgeColor','r');
PSNR4 = psnr(J_windows,J_refined);
imwrite(J_refined, 'J_refined.tif');

toc;




