
function morse(d)

  #Define constants
  De   = 11.2301 #eV
  beta =  2.627  #A^-1
  re   =  1.1282 #A
  r    = sqrt(d'd)

  #Solve for energy and force
  E = De * (1 - exp(-beta*r))^2
  F = 2 * De * (1 - 
end

function spring(d)

end


