rng shuffle;
clear; close all; clc;
% Create slider
%  sld = uicontrol('Style', 'slider',...
%      'Min',0.0,'Max',2.5,'Value',2.5,...
%      'Position', [400 20 120 20],...
%      'Callback', @surfzlim);
%% Config
% Position of robot and UAV
pos_2 = [2*cos(pi/6); 0; 0];
pos_3 = [0; 1; 0];
pos_1 = [0; -1; 0];
pos_robo = [pos_1 pos_2 pos_3];
mu_1 = [pos_1(1) pos_1(2)];
mu_2 = [pos_2(1) pos_2(2)];
mu_3 = [pos_3(1) pos_3(2)];


sigma_x1= [0.005 0; 0 0.005];
sigma_x2= [0.015 0; 0 0.015];
sigma_x3= [0.005 0; 0 0.005];

r_1 = mvnrnd(mu_1, sigma_x1, 1000);
r_2 = mvnrnd(mu_2, sigma_x2, 1000);
r_3 = mvnrnd(mu_3, sigma_x3, 1000);

% x1_mu=-0.03554+0.04221*x-0.001394*x^2 平均値と真値の差分(p.57)
% x2_mu=-0.04175-0.06191*x+0.0166*x^2
% x3_mu=-0.0008717-0.0007153*x+0.0003077*x^2

len_err = [2.14814214e-04  -2.98794337e-05  -2.07222678e-04];
height_err = 0.05;


% Prot result
font_size = 19;
markerSize = 30;
markerSize_uav = 10;

uav=[]

% Observe
% pos_uav = [1/cos(pi/6); 0; 2.0];　UAVの位置
% len_1u = norm(pos_1 - pos_uav);　normは幾何学的ベクトルの長さ(ベクトルの大きさ)
% len_2u = norm(pos_2 - pos_uav);
% len_3u = norm(pos_3 - pos_uav);
% len = [len_1u len_2u len_3u];

% FILE *fp; 
% 	char fname[] = "x1-trilateration.txt";
%     char readline[N] = {'\0'};
    
% type x1-trilateration.txt 

fileID_x1 = fopen('x1-trilateration.txt','r');
formatSpec = '%f';
x1 = fscanf(fileID_x1,formatSpec)
fclose(fileID_x1);

fileID_x2 = fopen('x2-trilateration.txt','r');
formatSpec = '%f';
x2 = fscanf(fileID_x2,formatSpec)
fclose(fileID_x2)

fileID_x3 = fopen('x3-trilateration.txt','r');
formatSpec = '%f';
x3 = fscanf(fileID_x3,formatSpec)
fclose(fileID_x3)


pos_uav=[cos(pi/6); 0; 2];
len = [norm(pos_1 - pos_uav) norm(pos_2 - pos_uav) norm(pos_3 - pos_uav)];
% len = [x1(v) x2(v) x3(v)];

std_devi = len_err(1) .* len.^2 + len_err(2) .* len + len_err(3);
measure_len = len + std_devi .* randn(1, 3);
measure_height = 2.0 + 2.0 * height_err .* randn(1, 1);


for j = 1:1:200
      samp_pos = [datasample(r_1, 1) 0;
                  datasample(r_2, 1) 0;
                  datasample(r_3, 1) 0]';

      len_std_devi = len_err(1) .* measure_len.^2 + len_err(2) .* measure_len + len_err(3);
      samp_len = measure_len + len_std_devi .* randn(1, 3);


      height_std_devi = measure_height * height_err;
      height_range = measure_height + height_std_devi .* [3; -3];

      [pos_u, signk3, k123, pos_b] = trilateration(samp_pos, samp_len);
      %[solvable, m, a, tangent] = has_solv(samp_pos, samp_len);
      
       plot3(pos_u(1), pos_u(2), pos_u(3),'.', 'Color','c');
       plot3(pos_u(1), pos_u(2), 0, '.', 'Color', '#808080');
       plot3(pos_u(1), 0.75, pos_u(3), '.', 'Color', '#808080');
       alpha(.5)
       hold on;
end

for v = 1:1:200

 axis equal;
% axis square;
grid on;


% len = [x1(v) x2(v) x3(v)];
% x1_mu=-0.03554+0.04221*(x1(v))-0.001394*(x1(v)^2)　平均値と真値の差分
% x2_mu=0.04175-0.06191*(x2(v))+0.0166*(x2(v)^2)
% x3_mu=0.0008717-0.0007153*(x3(v))+0.0003077*(x3(v)^2)
% len = [x1(v)-x1_mu x2(v)-x2_mu x3(v)-x3_mu];

 len = [x3(v)+0.05 x3(v) x3(v)-0.14];

% r1=len(1,1)
% r2=len(1,2)
% r3=len(1,3)

 [pos_u, signk3, k123, pos_b] = trilateration(pos_robo, len);

%  pos_u
%  k123

% pos_robo

% UGV plot
plot3(pos_robo(1,1),pos_robo(2,1),pos_robo(3,1),'.g', 'MarkerSize', markerSize)
text(pos_robo(1,1)+0.1,pos_robo(2,1)+0.1,pos_robo(3,1)+0.1,'P_1 Lens:x1','Color','g','FontSize',15,'FontAngle','italic')
hold on;

plot3(pos_robo(1,2),pos_robo(2,2),pos_robo(3,2),'.r', 'MarkerSize', markerSize)
text(pos_robo(1,2)+0.1,pos_robo(2,2)+0.1,pos_robo(3,2)+0.1,'P_2 Lens:x2','Color','r','FontSize',15,'FontAngle','italic')
hold on;

plot3(pos_robo(1,3),pos_robo(2,3),pos_robo(3,3),'.b', 'MarkerSize', markerSize)
text(pos_robo(1,3)+0.1,pos_robo(2,3)+0.1,pos_robo(3,3)+0.1,'P_3 Lens:x3','Color','blue','FontSize',15,'FontAngle','italic')
hold on;

% plot3(pos_robo(1, :), pos_robo(2, :), pos_robo(3, :), 'k.', 'MarkerSize', markerSize);
% 半球のplot

%plot3(pos_u(1), pos_u(2), pos_u(3), '+k','MarkerSize', markerSize_uav);

plot3(cos(pi/6), 0, 2, '.r','MarkerSize', markerSize);
text(1+0.15,0+0.15,1.5+0.15,'AprilTag','Color','black','FontSize',15,'FontAngle','italic')
% comet3(pos_u(1), pos_u(2), pos_u(3))
% text(pos_u(1)+0.1, pos_u(2)+0.1, pos_u(3)+0.1,'UAV','Color','black','FontSize',15)
hold on;


% uav_x=[pos_u(1);]
% uav_y=[pos_u(2);]
% uav_z=[pos_u(3);]
% 
% 
% fileID = fopen('uav_x.txt','a');
% % fprintf(fileID,'%6s %6s %6s\n','x','y','z');
% fprintf(fileID,'%12.8f \n',uav_x);
% fclose(fileID);
% 
% fileID = fopen('uav_y.txt','a');
% % fprintf(fileID,'%6s %6s %6s\n','x','y','z');
% fprintf(fileID,'%12.8f \n',uav_y);
% fclose(fileID);
% 
% fileID = fopen('uav_z.txt','a');
% % fprintf(fileID,'%6s %6s %6s\n','x','y','z');
% fprintf(fileID,'%12.8f \n',uav_z);
% fclose(fileID);
%  [x,y,z] = sphere(100);
%  s = surf(r1*x+pos_1(1,1),r1*y+pos_1(2,1),r1*z+pos_1(3,1),'FaceAlpha',0.15,'FaceColor','red','EdgeColor','none');
%  hold on;
%  
%  [x,y,z] = sphere(100);
%  s = surf(r2*x+pos_2(1,1),r2*y+pos_2(2,1),r2*z+pos_2(3,1),'FaceAlpha',0.15,'FaceColor','green','EdgeColor','none');
%  hold on;
%  
% [x,y,z] = sphere(100);
% s = surf(r3*x+pos_3(1,1),r3*y+pos_3(2,1),r3*z+pos_3(3,1),'FaceAlpha',0.15,'FaceColor','blue','EdgeColor','none');
% hold on;
% 
% 
% comet3(pos_u(1), pos_u(2), pos_u(3))
end

% zlim([0 5])
% 標準偏差プロット


%plot3(3, 0, 1.5, '.r','MarkerSize', markerSize);
%text(1+0.15,0+0.15,1.5+0.15,'Apriltag','Color','black','FontSize',15,'FontAngle','italic')
% comet3(pos_u(1), pos_u(2), pos_u(3))
% text(pos_u(1)+0.1, pos_u(2)+0.1, pos_u(3)+0.1,'UAV','Color','black','FontSize',15)
hold on;
  
F = getframe;
figure
imshow(F.cdata)
view(90,0)
xlabel('x [m]', 'FontSize', font_size);
ylabel('y [m]', 'FontSize', font_size);
zlabel('z [m]', 'FontSize', font_size);
uav
% 
% fileID = fopen('uav.txt','w');
% fprintf(fileID,'%6s %6s %6s\n','x','y','z');
% fprintf(fileID,'%12.8f %12.8f %12.8f\n',uav);
% fclose(fileID);
% 軌跡を取るならaxis equalでやる
%  axis equal;
% % axis square;
% grid on;
% hold off;


function surfzlim(source,event)
  font_size = 19;
  markerSize = 40;
  val = source.Value;
  pos_1 = [0; 0; 0];
  pos_2 = [2*cos(pi/6); 1; 0];
  pos_3 = [2*cos(pi/6); -1; 0];
  pos_robo = [pos_1 pos_2 pos_3];
  len = [val val val];

  [pos_u, signk3, k123, pos_b] = trilateration(pos_robo, len);
  plot3(pos_robo(1, :), pos_robo(2, :), pos_robo(3, :), 'k.', 'MarkerSize', markerSize);
  hold on;
  plot3(pos_u(1), pos_u(2), pos_u(3), ' ob');
  

  
  xlabel('x [m]', 'FontSize', font_size);
  ylabel('y [m]', 'FontSize', font_size);
  zlabel('z [m]', 'FontSize', font_size);
  
  signk3

  axis equal;
  grid on;
  hold off;
  
end

