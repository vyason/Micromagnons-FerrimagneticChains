# Declare Tensor operators for spin 3/2

using ITensors,ITensorMPS
import ITensors:op

# dimensionality
#--------------------------------------
ITensors.space(::SiteType"S=3/2") = 4

# Operator Sz
#--------------------------------------
ITensors.op(::OpName"Sz",::SiteType"S=3/2") =
[+3/2   0    0    0
  0    +1/2  0    0
  0     0   -1/2  0
  0     0    0   -3/2]

# Operator S+
#--------------------------------------
ITensors.op(::OpName"S+",::SiteType"S=3/2") =
[0  √3   0   0 
 0   0  √4   0
 0   0   0  √3
 0   0   0   0]

# Operator S-
#--------------------------------------
ITensors.op(::OpName"S-",::SiteType"S=3/2") =
[0   0   0   0
√3   0   0   0
 0  √4   0   0
 0   0  √3   0]

# z projection state vectors
#-----------------------------------------
ITensors.state(::StateName"+3/2",::SiteType"S=3/2") = [1, 0, 0, 0]          # |+3/2>
ITensors.state(::StateName"+1/2",::SiteType"S=3/2") = [0, 1, 0, 0]          # |+1/2>
ITensors.state(::StateName"-1/2",::SiteType"S=3/2") = [0, 0, 1, 0]          # |-1/2>
ITensors.state(::StateName"-3/2",::SiteType"S=3/2") = [0, 0, 0, 1]          # |-3/2>