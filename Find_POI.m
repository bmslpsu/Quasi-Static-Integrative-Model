function []=Find_POI()
%% this function reads the image saved in this folder and allows the user
%to acquire information from the graphs in the image

A = imread('master.img-011.jpg'); %reads the image
image(A)
% selecting the graphs of interest
disp('click on the green graph. Rotation angle')
[x1,y1] = getpts;       %rotation angle
disp('click on the blue graph. Stroke angle')
[x2,y2] = getpts; %stroke position angle
disp('click on the red graph. Deviation angle')
[x3,y3] = getpts;       %wing deviation angle

%% get the location of axis
%finds the axes range and translates from pixel to unit axis
disp('click on the upper limit and origin of y-axis')
[x4,y4] = getpts; %largest and smalles point on the y scale. Around 80 and -80
dely=80/(y4(1)-y4(2));

%% this gets the time axis
disp('click on the origin and desired end for the x-axis')
[x5,y5] = getpts;

t_axis=input('What is the length of the axis');
delx=t_axis/(x5(2)-x5(1));
%% transforms from image to graph coordinate
for i=1:length(x1)
    x1n(i)=x1(i)*delx;
    y1n(i)=(y1(i)-y4(2))*dely;
end
for i=1:length(x2)
    x2n(i)=x2(i)*delx;
    y2n(i)=(y2(i)-y4(2))*dely;
end
for i=1:length(x3)
    x3n(i)=x3(i)*delx;
    y3n(i)=(y3(i)-y4(2))*dely;
end

%%
figure
plot(x1n,y1n)
figure
plot(x2n,y2n)
figure
plot(x3n,y3n)
%% not always needed
%save('Angles.mat','x1n','y1n','x2n','y2n','x3n','y3n')
%% interpolation. required so all vectors have same length and points
[x11n, ia, ic] = unique(x1n);
y11n=y1n(ic:end-1);
%%
xmin=max([min(x11n) min(x2n) min(x3n)]);
xmax=min([max(x1n) max(x2n) max(x3n)]);
xx=xmin:0.001:xmax;
%% interpolation to make data smoother
yy1=spline(x11n,y11n,xx);
yy2=spline(x2n,y2n,xx);
yy3=spline(x3n,y3n,xx);

%% plots to show final graphs
plot(xx,yy1)
hold on
plot(xx,yy2)
hold on
plot(xx,yy3)
%% 
%save('AnglesInter.mat','xx','yy1','yy2','yy3')
