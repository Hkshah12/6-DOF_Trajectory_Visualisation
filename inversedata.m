initial_guess = [20 20 20 20 20 20]; % in degrees
initial_guess = initial_guess*0.0175; % degrees to radians

syms x y z real;

T_desired = [
    -0.5 -0.866 0 x;
    0.866 -0.5 0 y;
    0 0 1 z;
    0 0 0 1];

max_iterations = 200;

tolerance = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SAMPLING OF DESIRED TRAJECTORY%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nrows, ncols] = size(sample_points);
list_desired_joint_angles = zeros(nrows,6);
for i = 1:nrows
    row = sample_points(i,:);
    T_desired(1:3,4) = row
    sizeinitialguess = size(initial_guess)
    sizetdesired = size(T_desired)
    [output_joint_angles output_error]= newton_inverse_kinematics(T_desired, initial_guess, max_iterations, tolerance);
    output_joint_angles = output_joint_angles/0.0175; % convert result from radians to degrees
    list_desired_joint_angles(i,:) = output_joint_angles;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%VISUALIZATION%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %output_joint_angles = output_joint_angles*0.0175;
    %initial_guess = forward_kinematics(output_joint_angles)
end

function [jointangles, error] = newton_inverse_kinematics(T_desired, initial_guess, max_iterations, tolerance)
    % T_desired: Desired end-effector transformation matrix
    % initial_guess: Initial guess for joint angles
    % max_iterations: Maximum number of iterations
    % tolerance: Tolerance for convergence
    
    joint_angles = initial_guess;
    error = [inf inf inf];
    iteration = 0;
    
    while (abs(error(1)) > tolerance || abs(error(2)) > tolerance || abs(error(3)) > tolerance) iteration < max_iterations
        % Compute the forward kinematics for the current joint angles
        T_current = forward_kinematics(joint_angles)
        jointangles = double(joint_angles)
        % Compute the error between desired and current end-effector positions
        norm_error = norm(T_desired(1:3, 4) - T_current(1:3, 4));
        Tdesired = T_desired(1:3, 4);
        Tcurrent = T_current(1:3, 4);
        error = T_desired(1:3, 4) - T_current(1:3, 4);

        % Compute the Jacobian matrix
        J = compute_jacobian(joint_angles);

        % Compute the update for joint angles using Newton's method
        delta_theta = pinv(J) * (error);
        
        % Update joint angles
        joint_angles = double(joint_angles)
        s1 = size(joint_angles);
        s2 = size(delta_theta');
        joint_angles = joint_angles + delta_theta';
        updatedjointangles = double(joint_angles)
        iteration = iteration + 1;
    end
    disp('LOOP ENDED ONCE');
    if iteration >= max_iterations
        disp('Inverse kinematics did not converge within the maximum number of iterations.');
    end
end

function T = forward_kinematics(joint_angles)

   syms h1 h2 h3 h4 h5 h6 real;
   d1 = 89.2; %mm
   d4 = 109.3; %mm
   d5 = 94.75; %mm
   d6 = 82.5; %mm
   a2 = 425; %mm
   a3 = 392; %mm
   r11 = cos(h1)*cos(h2+h3+h4)*cos(h5)*cos(h6) + cos(h6)*sin(h1)*sin(h5) - cos(h1)*sin(h2+h3+h4)*sin(h6);
   r21 = cos(h2+h3+h4)*cos(h5)*cos(h6)*sin(h1) - cos(h1)*cos(h6)*sin(h5) - sin(h1)*sin(h2+h3+h4)*sin(h6);
   r31 = cos(h5)*cos(h6)*sin(h2+h3+h4) + cos(h2+h3+h4)*sin(h6);
   r12 = -cos(h1)*cos(h2+h3+h4)*cos(h5)*cos(h6) - sin(h1)*sin(h5)*sin(h6) - cos(h1)*cos(h6)*sin(h2+h3+h4);
   r22 = -cos(h2+h3+h4) + cos(h1)*sin(h5)*sin(h6) - cos(h6)*sin(h1)*sin(h2+h3+h4);
   r32 = -cos(h5)*sin(h2+h3+h4)*sin(h6) + cos(h2+h3+h4)*cos(h6);
   r13 = -cos(h1)*cos(h2+h3+h4)*sin(h5) + cos(h5)*sin(h1);
   r23 = -cos(h2+h3+h4)*sin(h1)*sin(h5) - cos(h1)*cos(h5);
   r33 = -sin(h2+h3+h4)*sin(h5);
   %p_x = r13*d6 + cos(h1)*(sin(h2+h3+h4)*d5 + cos(h2+h3)*a3 + cos(h2)*a2) + sin(h1)*d4;
   %p_y = r23*d6 + sin(h1)*(sin(h2+h3+h4)*d5 + cos(h2+h3)*a3 + cos(h2)*a2) - cos(h1)*d4;
   %p_z = r33*d6 - cos(h2+h3+h4)*d5 + sin(h2+h3)*a3 + sin(h2)*a2 + d1;
   p_x = d5*cos(h1)*sin(h2+h3+h4) + d4*sin(h1) - d6*cos(h1)*cos(h2+h3+h4) + a2*cos(h1)*cos(h2) + d6*cos(h5)*sin(h1) + a3*cos(h1)*cos(h2)*cos(h3) - a3*cos(h1)*sin(h2)*sin(h3);
   p_y = d5*sin(h1)*sin(h2+h3+h4) - d4*cos(h1) - d6*sin(h1)*cos(h2+h3+h4) - d6*cos(h1)*cos(h5) + a2*cos(h2)*sin(h1) + a3*cos(h2)*cos(h3)*sin(h1) - a3*sin(h1)*sin(h2)*sin(h3);
   p_z = d1 - d6*sin(h2+h3+h4)*sin(h5) + a3*sin(h2+h3) + a2*sin(h2) - d5*cos(h2+h3+h4);

   T = [
       r11 r12 r13 p_x;
       r21 r22 r23 p_y;
       r31 r32 r33 p_z;
       0 0 0 1];

   h1_value = joint_angles(1);
   h2_value = joint_angles(2);
   h3_value = joint_angles(3);
   h4_value = joint_angles(4);
   h5_value = joint_angles(5);
   h6_value = joint_angles(6);
   
   T_substituted = subs(T, h1, h1_value);
   T_substituted = subs(T_substituted, h2, h2_value);
   T_substituted = subs(T_substituted, h3, h3_value);
   T_substituted = subs(T_substituted, h4, h4_value);
   T_substituted = subs(T_substituted, h5, h5_value);
   T_substituted = subs(T_substituted, h6, h6_value);
   T_numeric = double(T_substituted);
   T = T_numeric;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function J = compute_jacobian(joint_angles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    syms h1 h2 h3 h4 h5 h6 real;
    d1 = 89.2; %mm
    d4 = 109.3; %mm
    d5 = 94.75; %mm
    d6 = 82.5; %mm
    a2 = 425; %mm
    a3 = 392; %mm
    
    %r13 = -cos(h1)*cos(h2+h3+h4)*sin(h5) + cos(h5)*sin(h1);
    %r23 = -cos(h2+h3+h4)*sin(h1)*sin(h5) - cos(h1)*cos(h5);
    %r33 = -sin(h2+h3+h4)*sin(h5);
    %p_x = r13*d6 + cos(h1)*(sin(h2+h3+h4)*d5 + cos(h2+h3)*a3 + cos(h2)*a2) + sin(h1)*d4;
    %p_y = r23*d6 + sin(h1)*(sin(h2+h3+h4)*d5 + cos(h2+h3)*a3 + cos(h2)*a2) - cos(h1)*d4;
    %p_z = r33*d6 - cos(h2+h3+h4)*d5 + sin(h2+h3)*a3 + sin(h2)*a2 + d1;

    %J_A = [
    %    0 sin(h1) sin(h1) sin(h1) cos(h1)*sin(h2+h3+h4) r13;
    %    0 -cos(h1) -cos(h1) -cos(h1) sin(h1)*sin(h2+h3+h4) r23;
    %    1 0 0 0 -cos(h2+h3+h4) r33];

    %J_L = [
    %    -p_y -cos(h1)*(p_z-d1) cos(h1)*(sin(h2+h3+h4)*sin(h5)*d6+cos(h2+h3+h4)*d5-sin(h2+h3)*a3) cos(h1)*(sin(h2+h3+h4)*sin(h5)*d6 + cos(h2+h3+h4)*d5) -d6*(sin(h1)*sin(h5) + cos(h1)*cos(h2+h3+h4)*cos(h5)) 0;
    %     p_x -sin(h1)*(p_z-d1) sin(h1)*(sin(h2+h3+h4)*sin(h5)*d6+cos(h2+h3+h4)*d5-sin(h2+h3)*a3) sin(h1)*(sin(h2+h3+h4)*sin(h5)*d6 + cos(h2+h3+h4)*d5) d6*(cos(h1)*sin(h5)-cos(h2+h3+h4)*cos(h5)*sin(h1)) 0;
    %     0 sin(h1)*p_y+cos(h1)*p_x -cos(h2+h3+h4)*sin(h5)*d6+sin(h2+h3+h4)*d5+cos(h2+h3)*a3 -cos(h2+h3+h4)*sin(h5)*d6+sin(h2+h3+h4)*d5 -cos(h5)*sin(h2+h3+h4)*d6 0];

    %J = [J_A;
    %    J_L];
    dp_x_1 = -d5*sin(h1)*sin(h2+h3+h4) + d4*cos(h1) + d6*sin(h1)*cos(h2+h3+h4) - a2*sin(h1)*cos(h2) + d6*cos(h5)*cos(h1) - a3*sin(h1)*cos(h2)*cos(h3) + a3*sin(h1)*sin(h2)*sin(h3);
    dp_x_2 = d5*cos(h1)*cos(h2+h3+h4) + d6*cos(h1)*sin(h2+h3+h4) - a2*cos(h1)*sin(h2) - a3*cos(h1)*sin(h2)*cos(h3) - a3*cos(h1)*cos(h2)*sin(h3);
    dp_x_3 = d5*cos(h1)*cos(h2+h3+h4) + d6*cos(h1)*sin(h2+h3+h4) - a3*cos(h1)*cos(h2)*sin(h3) - a3*cos(h1)*sin(h2)*cos(h3);
    dp_x_4 = d5*cos(h1)*cos(h2+h3+h4) + d6*cos(h1)*sin(h2+h3+h4);
    dp_x_5 = -d6*sin(h5)*sin(h1);
    dp_x_6 = 0;

    dp_y_1 = d5*cos(h1)*sin(h2+h3+h4) + d4*sin(h1) - d6*cos(h1)*cos(h2+h3+h4) + d6*sin(h1)*cos(h5) + a2*cos(h2)*cos(h1) + a3*cos(h2)*cos(h3)*cos(h1) - a3*cos(h1)*sin(h2)*sin(h3);
    dp_y_2 = d5*sin(h1)*cos(h2+h3+h4) + d6*sin(h1)*sin(h2+h3+h4) - a2*sin(h2)*sin(h1) - a3*sin(h2)*cos(h3)*sin(h1) - a3*sin(h1)*cos(h2)*sin(h3);
    dp_y_3 = d5*sin(h1)*cos(h2+h3+h4) + d6*sin(h1)*sin(h2+h3+h4) - a3*cos(h2)*sin(h3)*sin(h1) - a3*sin(h1)*sin(h2)*cos(h3);
    dp_y_4 = d5*sin(h1)*cos(h2+h3+h4) + d6*sin(h1)*sin(h2+h3+h4);
    dp_y_5 = d6*cos(h1)*sin(h5);
    dp_y_6 = 0;

    dp_z_1 = 0;
    dp_z_2 = -d6*cos(h2+h3+h4)*sin(h5) + a3*cos(h2+h3) + a2*cos(h2) + d5*sin(h2+h3+h4);
    dp_z_3 = -d6*cos(h2+h3+h4)*sin(h5) + a3*cos(h2+h3) + d5*sin(h2+h3+h4);;
    dp_z_4 = -d6*cos(h2+h3+h4)*sin(h5) + d5*sin(h2+h3+h4);
    dp_z_5 = -d6*sin(h2+h3+h4)*cos(h5);
    dp_z_6 = 0;

    J = [
        dp_x_1 dp_x_2 dp_x_3 dp_x_4 dp_x_5 dp_x_6;
        dp_y_1 dp_y_2 dp_y_3 dp_y_4 dp_y_5 dp_y_6;
        dp_z_1 dp_z_2 dp_z_3 dp_z_4 dp_z_5 dp_z_6];

    h1_value = joint_angles(1);
    h2_value = joint_angles(2);
    h3_value = joint_angles(3);
    h4_value = joint_angles(4);
    h5_value = joint_angles(5);
    h6_value = joint_angles(6);
    J_substituted = subs(J, h1, h1_value);
    J_substituted = subs(J_substituted, h2, h2_value);
    J_substituted = subs(J_substituted, h3, h3_value);
    J_substituted = subs(J_substituted, h4, h4_value);
    J_substituted = subs(J_substituted, h5, h5_value);
    J_substituted = subs(J_substituted, h6, h6_value);
    J_numeric = double(J_substituted);
    J = J_numeric;
    
end