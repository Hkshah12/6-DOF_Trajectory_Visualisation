
%Project2: ME5250
%Submitted by: Haard Shah
%Code: Visualization of Robot trajectory


joint_angles = d_list_angles * 0.01745; % Convert to radians

% Initialize the UR5 Denavit-Hartenberg (DH) Parameters
d1 = 162.5;  
a2 = 425;   
a3 = 392.2;   
d4 = 133.3; 
d5 = 99.7; 
d6 = 99.6;  

dhparams = [
    0, pi/2, L1, 0;  
    L2, 0, 0, 0;     
    L3, 0, 0, 0;     
    0, pi/2, L4, 0;  
    0, -pi/2, L5, 0; 
    0, 0, L6, 0;     
];

% Define the link lengths
links = [L1 L2 L3 L4 L5 L6];

% Create the figure and set the initial configuration
figure;
ax = gca; % Get current axis
hold on; 
plotUR5(joint_angles(1, :), dhparams, links, ax); % Display the initial configuration
axis equal; % Keep axis scaling consistent
%axis([-500 500 -500 500 -200 800])
view(3); % 3D view
title('UR5e Manipulator Robot Animation');
grid on;

% Set the animation speed (frames per second)
framesPerSecond = 5;
dt = 1 / framesPerSecond;

% Animate the robot motion using the joint angle list
%for i = 2:size(joint_angles, 1) % Start from the second set of joint angles
%    cla(ax); % Clear the plot to avoid overlapping
%    plotUR5(joint_angles(i, :), dhparams, links, ax); % Show new configuration
%    drawnow; % Refresh the plot
%    pause(dt); % Control frame rate
%end
% Animate the robot motion using the joint angle list
for i = 2:size(joint_angles, 1) % Start from the second set of joint angles
    cla(ax); % Clear the plot to avoid overlapping
    plotUR5(joint_angles(i, :), dhparams, links, ax); % Show new configuration
    % Update end-effector trajectory
    end_effector_trajectory(i, :) = getEndEffectorPosition(joint_angles(i, :), dhparams);
    % Plot trajectory
    plot3(end_effector_trajectory(1:i, 1), end_effector_trajectory(1:i, 2), end_effector_trajectory(1:i, 3), 'r:', 'LineWidth', 1);
    drawnow; % Refresh the plot
    pause(dt); % Control frame rate
end

% Function to plot the UR5 robot given the joint angles and DH parameters
function plotUR5(q, dhparams,links, ax)
    T = eye(4); % Initialize transformation matrix
    
    % Initialize end-effector position
    end_effector = zeros(3, size(q, 2));
    
    % Define colors for each link
link_colors = {'r', 'g', 'b', 'k', 'm', 'y'};

% Plot each link with its individual color
for i = 1:length(q)
    T = T * DHTransform(dhparams(i, :), q(i)); % Update transformation matrix
    
    % Plot link with individual color
    if i > 1
        plot3([T(1, 4) end_effector(1, i-1)], [T(2, 4) end_effector(2, i-1)], [T(3, 4) end_effector(3, i-1)], link_colors{i}, 'LineWidth', 3);
    end
    end_effector(:, i) = T(1:3, 4);
    hold on;
end

% Plot end-effector
plot3(T(1, 4), T(2, 4), T(3, 4), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

end

% Function to compute the Denavit-Hartenberg transformation matrix
function T = DHTransform(params, q)
    a = params(1);
    alpha = params(2);
    d = params(3);
    theta = params(4);
    
    T = [
        cos(q), -sin(q)*cos(alpha), sin(q)*sin(alpha), a*cos(q);
        sin(q), cos(q)*cos(alpha), -cos(q)*sin(alpha), a*sin(q);
        0, sin(alpha), cos(alpha), d;
        0, 0, 0, 1
    ];
end
% Function to compute the end-effector position given joint angles and DH parameters
function end_effector_position = getEndEffectorPosition(q, dhparams)
    T = eye(4);
    for i = 1:length(q)
        T = T * DHTransform(dhparams(i, :), q(i));
    end
    end_effector_position = T(1:3, 4)';
end
