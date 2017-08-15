% Analyze force and positions
%% Open tiff files and plot positions
FileInputName = 'Three dots.tiff';
info = imfinfo(FileInputName);
data=imread(FileInputName, 1, 'Info', info);
num_images = numel(info);       
dim = [info(1).Height info(1).Width num_images];
PixelSize = 106.7;  % PixelSize in nm
PointNo = 3;
%imshow(imadjust(data)); hold on;

% Import positions
Input = xlsread('Spot positions.xlsx');
Input(:,2:3) = Input(:,2:3)*PixelSize;
Tracks = cell(PointNo,1);
for i = 1:PointNo
    Tracks{i}=Input(Input(:,1)==(i-1),:); %Tracks{i}(:,2:3)=Tracks{i}(:,2:3)-Tracks{i}(1,2:3);
    plot(Tracks{i}(:,2),Tracks{i}(:,3)); hold on;
end

xline=[(min(Input(:,2))-((max(Input(:,2))-min(Input(:,2)))/10)) max(Input(:,2))];
dcm_obj = datacursormode;
%% Calculate Rotated Trace
Cursor = getCursorInfo(dcm_obj);
xPoints = [Cursor(1).Position(1) Cursor(2).Position(1)];
yPoints = [Cursor(1).Position(2) Cursor(2).Position(2)];
gradient = (yPoints(2)-yPoints(1))/(xPoints(2)-xPoints(1));
intercept = yPoints(1)-gradient*xPoints(1);
yline=gradient*xline+intercept;
plot(xline,yline,'red','LineWidth',2);
TracksRotated = Tracks;
Distance = Tracks;

figure;
% Find out rotated trace
for i = 1:PointNo
    x=Tracks{i}(:,2)-xline(1);
    y=Tracks{i}(:,3)-yline(1);
    radius=x.*x+y.*y;
    theta = atan(gradient);
    tanNew = tan(atan(y./x)-theta);
    TracksRotated{i}(:,2) = sqrt(radius./(1+tanNew.*tanNew));
    TracksRotated{i}(:,3) = tanNew.*TracksRotated{i}(:,2);
    Distance{i}(:,2:3)=TracksRotated{i}(:,2:3)-TracksRotated{i}(1,2:3);
    plot(TracksRotated{i}(:,2),TracksRotated{i}(:,3)); hold on;
    %plot(Distance{i}(:,2),Distance{i}(:,3)); hold on;
    dcm_obj = datacursormode;
end
%% Plot distances on graph for TracksRotated
figure;
MaxIndex = 0;
for i = 1:PointNo
    plot(TracksRotated{i}(:,4),TracksRotated{i}(:,2)); hold on;
    MaxIndex = max([MaxIndex max(TracksRotated{i}(:,4))]);
end
xlim([0 MaxIndex]);
% Add legend
TrackName = cell(PointNo,1);
for i = 1:PointNo
    TrackName{i} = ['Spot ' num2str(i)];
end
legend(TrackName{1},TrackName{2},TrackName{3});
dcm_obj = datacursormode;
%% Draw three lines based on average of six cursors. One line using two cursors
Cursor = getCursorInfo(dcm_obj);
for i = 1:PointNo
    % Find out the cursors for the particular track
    CursorData = zeros(2,1);
    CursorIndex = 1;
    for j = 1:PointNo*2
        if Cursor(j).Target.DisplayName==TrackName{i}
                CursorData(CursorIndex) = [Cursor(j).DataIndex];
                CursorIndex = CursorIndex + 1;
        end
    end
    CursorData = sort(CursorData);
    AverageLine = mean(TracksRotated{i}(CursorData(1):CursorData(2),2));
    % Plot line
    plot([TracksRotated{i}(1,4) TracksRotated{i}(end,4)],[AverageLine AverageLine],'-k'); hold on;
end

%  Save figure
delete(findall(gcf,'Type','hggroup'));
print('-djpeg','-r300','Three traces distances.jpg');
%% Plot Distances
figure;
for i = 1:PointNo
    plot(Distance{i}(:,2),Distance{i}(:,3)); hold on;
    dcm_obj = datacursormode;
end
%% Adjust center for Distance vector if necessary
i = 3;
Cursor = getCursorInfo(dcm_obj);
Distance{i}(:,2:3)=Distance{i}(:,2:3)-[Cursor(1).Position(1) Cursor(1).Position(2)];
figure;
for i = 1:PointNo
    plot(Distance{i}(:,2),Distance{i}(:,3)); hold on;
end
%% Plot overlapping distances on graph for Distance vector
figure;
for i = 1:3
    plot(Distance{i}(:,4),Distance{i}(:,2)); hold on;
end
% Add legend
TrackName = cell(PointNo,1);
for i = 1:PointNo
    TrackName{i} = ['Spot ' num2str(i)];
end
legend(TrackName{1},TrackName{2},TrackName{3});
print('-djpeg','-r300','Three traces distances.jpg');
%% Calculate and plot force
Force = Distance;
ForceSign = Distance;
L0 = 516;   % Contour length in nm
kT = 4.114; % Boltzmann constant time temperature in pN nm
P = 45;     % Persistence length in nm
figure;
Offset = 2; % Graph offset in pN
% Convert distance to force
for i=1:PointNo
    ForceSign{i}(:,2:3) = sign(Force{i}(:,2:3));
    Force{i}(:,2:3) = abs(Force{i}(:,2:3));
    [row,col] = find(Force{i}(:,2:3)>L0*0.90);
    Force{i}(row,col+1)=L0*0.90;

    Force{i}(:,2:3)=(1/4*(1-(Force{i}(:,2:3)/L0)).^(-2)-1/4+(Force{i}(:,2:3)/L0)-0.8*(Force{i}(:,2:3)/L0).^2.15)*kT/P;
    Force{i}(:,2:3)=Force{i}(:,2:3).*ForceSign{i}(:,2:3);
    plot(Force{i}(:,2),Force{i}(:,3)); hold on;
end

% Plot forces on graph
figure;
for i = 1:PointNo
    plot(Force{i}(:,4),Force{i}(:,2)+Offset*(i-1)); hold on;
end
TrackName = cell(PointNo,1);
MaxIndex = 0;
for i = 1:PointNo
    TrackName{i} = ['Spot ' num2str(i)];
    MaxIndex = max([MaxIndex max(TracksRotated{i}(:,4))]);
end
xlim([0 MaxIndex]);
legend(TrackName{1},TrackName{2},TrackName{3});
title('Force vs Time Graph');
xlabel('Time (frame)');
ylabel('Force (pN)');

%  Save figure
print('-djpeg','-r300','Forces Offset.jpg');
