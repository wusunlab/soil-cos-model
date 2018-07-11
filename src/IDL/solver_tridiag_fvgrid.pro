function solver_tridiag_FVgrid, C_0, alpha, beta, gamma, dt, Status_Check=status_check
; Tridiagonal linear equations solver for finite volume grid
; n - dimension of the linear system
; alpha, coefficients
; beta, coefficients
; gamma - source-sink term vector
; C - solution vector at the next time step
  n = n_elements(C_0);
  if (n_elements(alpha) ne n) or (n_elements(beta) ne n) or (n_elements(gamma) ne n) then $
    message, 'ERROR: dimensions of input vectors do not agree!'
    
    
  A_coef_mat = diag_matrix(alpha)
  ;A_coef_mat_inv = diag_matrix(1./alpha)
  
  beta_prime = beta
  beta_prime[0:n-2] = beta[1:n-1] 
  beta_prime[-1] = 0.
  B_coef_mat = - diag_matrix(beta) - diag_matrix(beta_prime)
  for loop_num = 1, n-1 do begin
    B_coef_mat[loop_num, loop_num-1] = beta[loop_num]
    B_coef_mat[loop_num-1, loop_num] = beta[loop_num]
  end

  S_coef_mat = gamma
  
  D_coef_mat = invert(2 * A_coef_mat - dt * B_coef_mat) 
    
  C = D_coef_mat ## ( (2 * A_coef_mat + dt * B_coef_mat) ## C_0 + 2 * dt * S_coef_mat )
  
  return, C
end