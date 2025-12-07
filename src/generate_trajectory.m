% Define the size of the square (in millimeters)
square_size = 100; % Assuming a 100x100 mm square

% Define the interval (in millimeters)
interval = 1; % 1mm interval

% Compute the number of points per side
num_points_per_side = floor(square_size / interval);

% Initialize an array to store the sampled points
points = [];

% Sample points along the top side
points = [points; [0:interval:square_size]', ones(num_points_per_side + 1, 1) * square_size];

% Sample points along the right side (excluding the top-right corner)
points = [points; ones(num_points_per_side, 1) * square_size, (square_size - interval:-interval:0)'];

% Sample points along the bottom side (excluding the bottom-right corner)
points = [points; (square_size - interval:-interval:0)', zeros(num_points_per_side, 1)];

% Sample points along the left side (excluding the bottom-left corner)
points = [points; zeros(num_points_per_side, 1), (0:interval:square_size - interval)', ];

points = points + 50;
points(:,3) = 0;

% Display the sampled points
scatter(points(:,1), points(:,2), 'filled');
axis equal; % Equal aspect ratio
xlim([-1 square_size+1]);
ylim([-1 square_size+1]);
xlabel('X (mm)');
ylabel('Y (mm)');
title('Sampled Points on the Perimeter of a Planar Square');
