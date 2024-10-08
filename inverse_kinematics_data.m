
%Project2: ME5250
%Submitted by: Haard Shah
%Code: Calculating the Inverse Kinematics based on sampled trajectory

% guess_in: Initial guess for joint angles
%This is for the computation to begin as without this Newton Raphson cannot
%proceed without inital guess.
guess_in = [15 20 25 10 15 12]; % (units: in degrees)
guess_in = guess_in*0.0175; % (Conversion degrees to radians = pi/180 = 0.01745 ~ 0.0175)

syms a b c real;

%Desired Transformation Matrix
Desired_tr = [
    -0.5 -0.866 0 a;
    0.866 -0.5 0 b;
    0 0 1 c;
    0 0 0 1];

max_iterations = 300;
threshold = 1; %Upper limit of convergence

[nrows, ncols] = size(sample_points);
d_list_angles = zeros(nrows,6);
for i = 1:nrows
    row = sample_points(i,:);
    Desired_tr(1:3,4) = row;
    [out_joint, output_error]= newton_inverse_kinematics(Desired_tr, guess_in, max_iterations, threshold);
    out_joint = out_joint/0.0175; 
    d_list_angles(i,:) = out_joint;
    out_joint = out_joint*0.01744;
    guess_in = out_joint;

end

function [jointangles, error] = newton_inverse_kinematics(Desired_tr, guess_in, max_iterations, threshold)
      
    joint_angles = guess_in;
    error = [inf inf inf];
    iteration = 0;
    
    while (abs(error(1)) > threshold || abs(error(2)) > threshold || abs(error(3)) > threshold) && iteration < max_iterations
        % Compute the forward kinematics for the current joint angles
        T_current = forward_kinematics(joint_angles);
        jointangles = double(joint_angles);

        % Compute the error between desired and current end-effector positions
                
        error = Desired_tr(1:3, 4) - T_current(1:3, 4);

        % Compute the Jacobian matrix
        J = compute_jacobian(joint_angles);

        % Compute the update for joint angles using Newton's method
        delta_theta = pinv(J) * (error);
        
        % Update joint angles  
        joint_angles = joint_angles + delta_theta';
        
        iteration = iteration + 1;
    end
    disp('Loop is completed');
    if iteration >= max_iterations
        disp('Could not converge within the specified threshold.');
    end
end

function T = forward_kinematics(joint_angles)

   syms h1 h2 h3 h4 h5 h6 real;
   d1 = 162.5; 
   a2 = 425; 
   a3 = 392.2;
   d4 = 133.3; 
   d5 = 99.7; 
   d6 = 99.6; 
 

   e11 = cos(h1)*cos(h2+h3+h4)*cos(h5)*cos(h6) + cos(h6)*sin(h1)*sin(h5) - cos(h1)*sin(h2+h3+h4)*sin(h6);
   e12 = -cos(h1)*cos(h2+h3+h4)*cos(h5)*cos(h6) - sin(h1)*sin(h5)*sin(h6) - cos(h1)*cos(h6)*sin(h2+h3+h4);
   e13 = -cos(h1)*cos(h2+h3+h4)*sin(h5) + cos(h5)*sin(h1);
   e21 = cos(h2+h3+h4)*cos(h5)*cos(h6)*sin(h1) - cos(h1)*cos(h6)*sin(h5) - sin(h1)*sin(h2+h3+h4)*sin(h6);
   e22 = -cos(h2+h3+h4) + cos(h1)*sin(h5)*sin(h6) - cos(h6)*sin(h1)*sin(h2+h3+h4);
   e23 = -cos(h2+h3+h4)*sin(h1)*sin(h5) - cos(h1)*cos(h5);
   e31 = cos(h5)*cos(h6)*sin(h2+h3+h4) + cos(h2+h3+h4)*sin(h6); 
   e32 = -cos(h5)*sin(h2+h3+h4)*sin(h6) + cos(h2+h3+h4)*cos(h6);
   e33 = -sin(h2+h3+h4)*sin(h5);
 
   posx = d5*cos(h1)*sin(h2+h3+h4) + d4*sin(h1) - d6*cos(h1)*cos(h2+h3+h4) + a2*cos(h1)*cos(h2) + d6*cos(h5)*sin(h1) + a3*cos(h1)*cos(h2)*cos(h3) - a3*cos(h1)*sin(h2)*sin(h3);
   posy = d5*sin(h1)*sin(h2+h3+h4) - d4*cos(h1) - d6*sin(h1)*cos(h2+h3+h4) - d6*cos(h1)*cos(h5) + a2*cos(h2)*sin(h1) + a3*cos(h2)*cos(h3)*sin(h1) - a3*sin(h1)*sin(h2)*sin(h3);
   posz = d1 - d6*sin(h2+h3+h4)*sin(h5) + a3*sin(h2+h3) + a2*sin(h2) - d5*cos(h2+h3+h4);

   T = [
       e11 e12 e13 posx;
       e21 e22 e23 posy;
       e31 e32 e33 posz;
       0 0 0 1];

   h1_val = joint_angles(1);
   h2_val = joint_angles(2);
   h3_val = joint_angles(3);
   h4_val = joint_angles(4);
   h5_val = joint_angles(5);
   h6_val = joint_angles(6);
   
   T_sub = subs(T, h1, h1_val);
   T_sub = subs(T_sub, h2, h2_val);
   T_sub = subs(T_sub, h3, h3_val);
   T_sub = subs(T_sub, h4, h4_val);
   T_sub = subs(T_sub, h5, h5_val);
   T_sub = subs(T_sub, h6, h6_val);
   T_numeric = double(T_sub);
   T = T_numeric;
    
end

function J = compute_jacobian(joint_angles)

    syms h1 h2 h3 h4 h5 h6 real;
    d1 = 162.5; 
    d4 = 133.3; 
    d5 = 99.7; 
    d6 = 99.6; 
    a2 = 425; 
    a3 = 392.2; 
        
    dposx_1 = -d5*sin(h1)*sin(h2+h3+h4) + d4*cos(h1) + d6*sin(h1)*cos(h2+h3+h4) - a2*sin(h1)*cos(h2) + d6*cos(h5)*cos(h1) - a3*sin(h1)*cos(h2)*cos(h3) + a3*sin(h1)*sin(h2)*sin(h3);
    dposy_1 = d5*cos(h1)*sin(h2+h3+h4) + d4*sin(h1) - d6*cos(h1)*cos(h2+h3+h4) + d6*sin(h1)*cos(h5) + a2*cos(h2)*cos(h1) + a3*cos(h2)*cos(h3)*cos(h1) - a3*cos(h1)*sin(h2)*sin(h3);
    dposz_1 = 0;
    
    dposx_2 = d5*cos(h1)*cos(h2+h3+h4) + d6*cos(h1)*sin(h2+h3+h4) - a2*cos(h1)*sin(h2) - a3*cos(h1)*sin(h2)*cos(h3) - a3*cos(h1)*cos(h2)*sin(h3);
    dposy_2 = d5*sin(h1)*cos(h2+h3+h4) + d6*sin(h1)*sin(h2+h3+h4) - a2*sin(h2)*sin(h1) - a3*sin(h2)*cos(h3)*sin(h1) - a3*sin(h1)*cos(h2)*sin(h3);
    dposz_2 = -d6*cos(h2+h3+h4)*sin(h5) + a3*cos(h2+h3) + a2*cos(h2) + d5*sin(h2+h3+h4);

    dposx_3 = d5*cos(h1)*cos(h2+h3+h4) + d6*cos(h1)*sin(h2+h3+h4) - a3*cos(h1)*cos(h2)*sin(h3) - a3*cos(h1)*sin(h2)*cos(h3);
    dposy_3 = d5*sin(h1)*cos(h2+h3+h4) + d6*sin(h1)*sin(h2+h3+h4) - a3*cos(h2)*sin(h3)*sin(h1) - a3*sin(h1)*sin(h2)*cos(h3);
    dposz_3 = -d6*cos(h2+h3+h4)*sin(h5) + a3*cos(h2+h3) + d5*sin(h2+h3+h4);

    dposx_4 = d5*cos(h1)*cos(h2+h3+h4) + d6*cos(h1)*sin(h2+h3+h4);
    dposy_4 = d5*sin(h1)*cos(h2+h3+h4) + d6*sin(h1)*sin(h2+h3+h4);
    dposz_4 = -d6*cos(h2+h3+h4)*sin(h5) + d5*sin(h2+h3+h4);

    dposx_5 = -d6*sin(h5)*sin(h1);
    dposy_5 = d6*cos(h1)*sin(h5);
    dposz_5 = -d6*sin(h2+h3+h4)*cos(h5);

    dposx_6 = 0;    
    dposy_6 = 0;
    dposz_6 = 0;

    J = [
        dposx_1 dposx_2 dposx_3 dposx_4 dposx_5 dposx_6;
        dposy_1 dposy_2 dposy_3 dposy_4 dposy_5 dposy_6;
        dposz_1 dposz_2 dposz_3 dposz_4 dposz_5 dposz_6];

    h1_val = joint_angles(1);
    h2_val = joint_angles(2);
    h3_val = joint_angles(3);
    h4_val = joint_angles(4);
    h5_val = joint_angles(5);
    h6_val = joint_angles(6);

    J_sub = subs(J, h1, h1_val);
    J_sub = subs(J_sub, h2, h2_val);
    J_sub = subs(J_sub, h3, h3_val);
    J_sub = subs(J_sub, h4, h4_val);
    J_sub = subs(J_sub, h5, h5_val);
    J_sub = subs(J_sub, h6, h6_val);
    J_numeric = double(J_sub);
    J = J_numeric;
    
end