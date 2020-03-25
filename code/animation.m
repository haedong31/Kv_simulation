clc
close all
clear variables

%% AP & states associated with Ito
[t, S, A, C] = Rasmusson_AP(50);

plot(t, S(:,1))
% plot(t, A(:,69))
% plot(t, A(:,70))

hold on
p = plot(t(1), S(1,1), 'r*', 'MarkerFaceColor', 'red');
hold off

for k = 2:length(t)
    p.XData = t(k);
    p.YData = S(k,1);
    drawnow
end


%%
clc
close all
clear variables

I = imread('./graph.jpg');
gray_I = rgb2gray(I);

figure(1)
imshow(gray_I, [])
% d = imdistline;

[centers, radii] = imfindcircles(gray_I,[100,170],'ObjectPolarity','dark','Sensitivity',0.95);
