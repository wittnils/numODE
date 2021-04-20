function sheet01
  % t is a vector with 100 linearly spaced entries between 0 and 300 
  t = linspace (0 ,300 ,20);
  figure
  % hold on -> plotting commands wont be executed immediately
  hold on
  for i = 1:1
    % lsode: function to solve a system of ODEs. 
    % "lotka" is the name of system of functions on the RHS of dx/dt = f(x,t) 
    % Signature: lsode("NAME", x_0, t, t_crit), where x_0 is the initival value for the IVP
    % t is a vector that contains the time values, at which a solution is calculated
    x = lsode ("lotka" , [ 0.5 * i , 1 ] , t ) ;
    % for every time-step there is a row and in each row there are two entries: one for each equation containing the 
    % numerical solution of our system of ODEs at the given time-step
    % plot first and second column of x 
    plot (x (: ,1) , x ( : , 2 ) ) ;
  end
  x
  xlabel ( " x_1 " ) ;
  ylabel ( " x_2 " ) ;
endfunction

function xdot = lotka ( x , t )
  % Parameters for System of ODEs
  a = 1;
  b = 1;
  c = 1;
  d = 1;
  xdot = zeros(2,1) ;
  % Lotka-Volterra System of ODEs 
  xdot ( 1 ) = a*x(1) -b*x (1) * x (2) ;
  xdot ( 2 ) = c *x (1) * x(2) -d*x (2) ;
end