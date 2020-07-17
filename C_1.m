
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Target Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Anonymous function specification of F:R^2 -> R^2 
%
%            F1(v(1),v(2))       F2(v(1),v(2)) 
%                   |                  |     
F  = @(v)[20/(pi*(1+(20*v(2)-100)^2)), 1]';
%     			 
% The function returns a 2 x 1 column vector (note the ') 
%
%
tInitial    = 0.0;                         % Initial time
tFinal      = 10.0;                        % Final time
yInitial    = [(1/pi)*atan(-100)+0.5,0]';  % Initial value of y (2 x 1 column vector) 
E_glob_tol  = 1*10^-5;                     % Global error tolerance 

y           = zeros(2,1);                  % Arrays to hold solution values
LTE_est     = zeros(1,1);
error       = zeros(1,1);
t           = zeros(1,1);
i           = 1;

%
% ////////////////////////////////////////////////////////////
%    Computing the solution using adaptive timestep  method 
% ////////////////////////////////////////////////////////////
%

t(1)    = tInitial;
y(:,1)  = yInitial;  
Glob_err_est(1) = 0;
error(1) = 0;

% pick initial h so that |y_1 - y_0| is small
h = 0.1;

tic

while t(i) < tFinal
    % create temporary h 
    h_temp = h;
    
    % calculating next value with h_temp 
    t_hold       = t(i) + h_temp;
    y_hold       = y(:,i) + (h_temp/2)*F(y(:,i)) + (h_temp/2)*F(y(:,i)+h_temp*F(y(:,i)));
        
    % calculate next value with h_temp/2
    y_intermediate  = y(:,i) + (h_temp/4)*F(y(:,i)) + (h_temp/4)*F(y(:,i)+(h_temp/2)*F(y(:,i)));
    y_half_step = y_intermediate + (h_temp/4)*F(y_intermediate) + (h_temp/4)*F(y_intermediate + (h_temp/2)*F(y_intermediate));
    
    % calculate Err(h_temp)
    Err_h_temp = 2*abs(y_hold(1,1) - y_half_step(1, 1))/3;
    
    % calculate estimate to global error/time interval
    E_glob_est = Err_h_temp/h_temp;
    
    % continue reducing h_temp until restriction is satisfied 
    while E_glob_est > (E_glob_tol/(tFinal - tInitial))
        
        h_temp_hold = h_temp;
        h_temp = 0.95 * h_temp * sqrt((E_glob_tol/10)*(h_temp/Err_h_temp));
        
        % calculating next value with h_temp 
        t_hold       = t(i) + h_temp;
        
        % calculate new Err(h_temp)
        Err_h_temp = Err_h_temp*(h_temp/h_temp_hold)^3;
    
        % calculate estimate to global error/time interval
        E_glob_est = Err_h_temp/h_temp;
        
    end
    
    % if timestep is too small we increase it 
    while E_glob_est < 0.8*(E_glob_tol/(tFinal - tInitial))
        
        h_temp_hold = h_temp;
        h_temp = 1.05 * h_temp * sqrt((E_glob_tol/10)*(h_temp/Err_h_temp));
        
        % calculating next value with h_temp 
        t_hold       = t(i) + h_temp;
        
        % calculate new Err(h_temp)
        Err_h_temp = Err_h_temp*(h_temp/h_temp_hold)^3;
    
        % calculate estimate to global error/time interval
        E_glob_est = Err_h_temp/h_temp;
    end
    
    
    % finalise h once restriction is satisfied  
    h = h_temp;
    
    % calculate next t and y values with finalised h 
    t(i+1) = t(i) + h;
    y(:,i+1) = y_hold;
    
    % calculate exact global error and global error estimate
    error(i+1)   = abs(y(1,i+1) - exact_sol(t(i+1)));
    Glob_err_est(i+1) = E_glob_est * 10;
    
    % update i to i+1
    i = i+1;
   
end

toc

%    Set plot limits 
yMin =    0.0;
yMax =    2.0;

% plot solution, global error and global error estimate
figure(1)
plot(t,y(1,:));
axis([0,tFinal,yMin,yMax]);
 
figure(2);
plot(t, error)
figure(3);
plot(t, Glob_err_est);


% function for exact solution

function x = exact_sol(t)
    x = (1/pi)*atan(20*(t-5)) + 0.5;
end




