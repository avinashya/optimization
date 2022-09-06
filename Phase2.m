clc;
clear all;
que=2;
x0 = [-10;-6;8] 
optimumvalue = project_phase_2(x0,que)

function answer = project_phase_2(x0,que)
 % x0 = Initial Approximation
 k = 1 ; % Iteration Counter
 e = 0.001 ; % Accuracy
 fun = [] ;
while true
 gradient_0 = grad(x0,que) ;  % Gradient at initial Point
if magnitude(gradient_0) < e % Termination condition 1 
    answer = x0 ;
    
    break
end
 fun(k) = Function(x0,que)
hessian = hess(x0,que) ; % Hessian matrix initial Point
direction = inv(hessian)*gradient_0 ;% Direction for unidirectional Search 
Alpha = project_1(x0,direction,que) ; % Value of Alpha by unidirectional Search using Project Phase 1
x1 = new_point(x0,Alpha,direction) ;% New point after substituting the value
gradient_1 = grad(x1,que) ;
if (gradient_0*gradient_1' < e ) % Termination Condition 2 
    answer = x1 ;
    break
end
if ((magnitude(x1 - x0))/(magnitude(x0))) < e % Termination COndition 3
answer = x1;
   break
end
 
k = k + 1;
x0 = x1;
end
iteration = [1:k]
plot(iteration,fun)
end
% ------- Multivariable Function ------------
function fun_val = Function(x,que)
fun_val=0;
if (que==1)
    n=5;
for i=1:n
    fun_val=fun_val+(i*x(i))^2;
end
elseif (que==2)
    n=3;
    for i=1:(n-1)
    fun_val= fun_val + 100*(x(i+1)-(x(i))^2)^2 + (x(i)-1)^2;
    end
elseif (que==3)
    n=4;
    fun_val= (x(1)-1)^2;
    for i=2:n
    fun_val= fun_val + i*(2*(x(i))^2 - x(i-1))^2 ;
    end
elseif (que==4)
    n=6;
    fun_val=(x(1)-1)^2;
    for i=2:n
      fun_val=fun_val + (x(i)-1)^2 - x(i)*x(i-1);
    end
  
elseif (que==5)
    n=2;
    a=0;
    b=0;
    c=0;
    for i=1:n
        a= a + (x(i))^2;
        b= b + 0.5*i*x(i);
        c= c + 0.5*i*x(i);
    end
    fun_val= a + b^2 + c^4;
 
else
    fun_val=0;
end
 
end

%-----Function for gradient ----------------
function gradient = grad(x,que) 
gradient = zeros(length(x),1);
h = 0.001;
for i = 1:length(x)
    y = x;
 y(i) = y(i)+h;
 a = Function(y,que);
 y(i) = y(i)-2*h;
 b = Function(y,que);
 gradient(i) = (a - b)/(2*h);
end
end

%---Fuction creation for single variable optimization for gradient---
function fun_val = Objective_Fun(y,x0,direction,que)
l = length(x0);
for i = 1:l
x0(i) = x0(i) - y*direction(i);
end
fun_val = Function(x0,que);
end
%-------------Function for magnitude of a vector--------------
function m = magnitude(gradient)
magnitude_squeuare = gradient.*gradient;
magnitude_squeuare_sum = sum(magnitude_squeuare);
m = sqrt(magnitude_squeuare_sum);
end

%------------Function for Hessian matrix -------------------
function hessian = hess(x, que)
l = length(x) ;
hessian = zeros(l,l);
h = 0.001 ;
for i = 1:l
for j = 1:l
if i == j 
 y = x;
 y(i) = y(i)+h ;
 a = Function(y,que);
 y(i) = y(i)-2*h ;
 b = Function(y,que);
 c = Function(x,que);
 hessian(i,j) =  (a+b-2*c)/(h^2);
else
 a = x ;
 b = x ;
 c = x ;
 d = x ;
 a(i) = a(i) + h;
 a(j) = a(j) + h;
 first_term = Function(a,que);
 b(i) = b(i)+ h ;
 b(j) = b(j)- h ;
 second_term = Function(b,que);
 c(i) = c(i) - h ;
 c(j) = c(j) + h ;
 third_term = Function(c,que);
 d(i) = d(i) - h ;
 d(j) = d(j) - h ;
 forth_term = Function(d,que);
hessian(i,j)= (first_term - second_term - third_term + forth_term)/(4*h^2);
end
end
end
end

%------Function to find new vector after finding the value of Alpha-------
function new_points = new_point(x,Alpha,direction)
l = length(x) ;
for i = 1:l
    x(i) = x(i) - Alpha*direction(i) ;
end
new_points = x ;
end
%------------Function of single variable optimization----------
%------------Same as Project Phase 1 --------------------------
function Alpha = project_1(x,direction,que)
% couple code for phase 1
% Bounding phase
x0 = 0.6; % initial guess
delta = 0.5; % small increment
iter = 0;
xd1 = x0-abs(delta);
xd3 = x0+abs(delta);
fd1= Objective_Fun(xd1,x,direction,que); fd3= Objective_Fun(xd3,x,direction,que); f1 = Objective_Fun(x0,x,direction,que);
if (fd1 >= f1) && (f1 >= fd3) % if ok then delta +ve otherwise -ve
    delta = + abs(delta);
else
    delta = - abs(delta);
end
xi = x0;
f3 = 1;
f2 =2;
while f3 < f2 %termination condition
    f2 = Objective_Fun(xi,x,direction,que); % f(x(n)) function value at x new xn=xk
    x1 = xi + (2^iter)*delta; % x new
    f3 = Objective_Fun(x1,x,direction,que); % f(x(k+1)) function value at x(k+1)
	xh=xi;
    xi = x1;
    iter = iter +1; % next iteration
end
% Newton Raphson Method

delta = 0.01;
e = 1;
iter = 1;
xi = x1;
while e > 0.0001
    xd1 = xi + delta;
    xd2 = xi - delta;
    dfn =(Objective_Fun(xd1,x,direction,que)- Objective_Fun(xd2,x,direction,que))/(2*delta);
    sfn =(Objective_Fun(xd1,x,direction,que)- 2*Objective_Fun(xi,x,direction,que) + Objective_Fun(xd2,x,direction,que))/(delta)^2;
    xn = xi -(dfn)/(sfn);
    xd1 = xn + delta;
    xd2 = xn - delta;
    e  = abs((Objective_Fun(xd1,x,direction,que)- Objective_Fun(xd2,x,direction,que))/(2*delta));
    xi =xn;
    
end
Alpha = xn;
end  
   
