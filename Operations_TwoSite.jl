# Useful functions for two-site Tensor operations

using ITensors,ITensorMPS
import ITensors:op


#   entanglement negativity between two sites
# wf = input wavefunction in MPS form
# locs = [A,B] = end locations of the two sites in the chain
#-----------------------------------
function calc_Neg_TwoSite(wf::MPS,locs)
    psi = copy(wf)
    A,B = locs

    rho = TwoSiteDM(psi,[A,B])
    old = inds(rho)

    new = [old[1],old[2],old[4],old[3]]
    rho_PT = swapinds(rho,old,new)
    egn_val = calc_Eigenvalues(rho_PT)

    Neg =  abs(sum( [ lam for lam in egn_val if lam<0 ]))
    #-----------------------------------
    return Neg
end


#   mutual information between two sites
# wf = input wavefunction in MPS form
# locs = [A,B] = end locations of the two sites in the chain
#-----------------------------------
function calc_MutInf_TwoSite(wf::MPS,locs)
    psi = copy(wf)
    A,B = locs
    
    S_A  = calc_SvN( TwoSiteDM(psi,[A,A]) )
    S_B  = calc_SvN( TwoSiteDM(psi,[B,B]) )
    S_AB = calc_SvN( TwoSiteDM(psi,[A,B]) )

    MutInf = S_A + S_B - S_AB
    #-----------------------------------
    return MutInf
end


# returns a two sites reduced density matrix
# wf = input wavefunction in MPS form
# locs = [A,B] = end locations of the two sites in the chain
#--------------------------------------
function TwoSiteDM(wf::MPS,locs)
    ket = copy(wf)
    sites = siteinds(ket)

    A,B = locs
    orthogonalize!(ket,A)
    bra = prime(dag(ket),linkinds(ket))
    
    rho = prime(ket[A],linkinds(ket,A-1)) * prime(bra[A],sites[A])
    for k in A+1:B-1
        rho *= ket[k]*bra[k]
    end
    rho *= prime(ket[B],linkinds(ket,B))  * prime(bra[B],sites[B])  
    
    #-----------------------------------
    return rho
end