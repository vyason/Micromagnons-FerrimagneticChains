# Useful functions for four-site Tensor operations

using ITensors,ITensorMPS
import ITensors:op


# four point RDM of sites [A,A+1] [B,B+1]
# wf = input wavefunction in MPS form
# locs = [A,B] = locations of the sites in the chain
#------------------------------------------
function FourSiteDM(wf::MPS,locs)
    psi = copy(wf)
    A,C = locs
    sep = C-A-2

    orthogonalize!(psi,A)
    ket = psi[A]
    for k in A+1:C+1
        ket *= psi[k]
    end
    rho = prime(ket,"Site") * dag(ket)

    if sep>0
        inds_list = inds(rho)
        set1 = collect(3:C-A)
        set2 = set1[end] .+ set1 .+ 2
        for idx in 1:length(set1)
            rho = rho * delta(inds_list[set1[idx]],inds_list[set2[idx]])
        end
    end
    #-------------------
    return rho
end



#   4-partite entanglement details between four sites
# wf = input wavefunction in MPS form
# locs = [A,B] = locations of the sites in the chain
#-----------------------------------
function calc_LogNeg_FourSite(wf::MPS,locs)   

    #-----------------------------------
    function get_LogNeg_OnePartition(rho,old,new)
        rho_PT = swapinds(rho,old,new)
        egn_val = calc_Eigenvalues(rho_PT)
        Neg =  abs( sum( [lam for lam in egn_val if lam<0] ) )
        return log(1+2*Neg)
    end
    #-----------------------------------
    function exchange_indices(old,set1,set2)
        new = copy(old)
        for idx in 1:length(set1)
            new[set1[idx]],new[set2[idx]] = new[set2[idx]],new[set1[idx]]
        end
        return new
    end
    #-----------------------------------

    psi = copy(wf)
    A,C = locs

    rho = FourSiteDM(psi,[A,C])
    old = [inds(rho)...]

    #   ABCD = 1 2 3 4
    #   ABCD = 5 6 7 8
    new_list    = Vector{Any}(undef,7)
    new_list[1] = exchange_indices(old,[1],[5])         #A
    new_list[2] = exchange_indices(old,[2],[6])         #B
    new_list[3] = exchange_indices(old,[3],[7])         #C
    new_list[4] = exchange_indices(old,[4],[8])         #D
    new_list[5] = exchange_indices(old,[1,2],[5,6])     #AB
    new_list[6] = exchange_indices(old,[1,3],[5,7])     #AC
    new_list[7] = exchange_indices(old,[1,4],[5,8])     #AD

    LogNeg_Partitions = [ get_LogNeg_OnePartition(rho,old,new) for new in new_list ]
    LogNeg = prod(LogNeg_Partitions)^(1.0/7.0)

    #-----------------------------------
    return LogNeg
end