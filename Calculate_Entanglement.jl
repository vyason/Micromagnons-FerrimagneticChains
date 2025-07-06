using ITensors,ITensorMPS
import ITensors:op

#   calc_Eigenvalues of density matrix
#-----------------------------------
function calc_Eigenvalues(density_matrix)
    egn_val,_ = eigen(density_matrix,ishermitian=true)
    #-----------------------------------
    return diag(array(egn_val))
end


#   eigenvalue based generic properties of a density matrix
#-----------------------------------
function calc_Norm(density_matrix)
    rho = copy(density_matrix)

    egn_val = calc_Eigenvalues(rho)
    Norm = sum(egn_val)
    Purity = sum(egn_val.^2)

    #-----------------------------------
    return [Norm,Purity]
end


#   eigenvalue based generic properties of a density matrix
#-----------------------------------
function calc_SvN(density_matrix)
    rho = copy(density_matrix)

    egn_val = calc_Eigenvalues(rho)
    SvN = sum( [ - lam*log(lam) for lam in egn_val if lam > 0 ] )

    #-----------------------------------
    return SvN
end