% Computing the first order derivative for Newton-Raphson Function
function out = F_prime(f,x,dx)
  out = (F(f,x+dx) - F(f,x))/(dx)
endfunction
