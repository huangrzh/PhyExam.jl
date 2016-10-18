# --------------module PhyExam-----------
# huangrzh@icloud.com
#---------------------------------------------------



module PhyExam

# package code goes here
export ED_1D_hessenberg

type Ops
  eye::Array{Float64,2}
  sp::Array{Float64,2}
  sz::Array{Float64,2}
  heff::Array{Float64,2}
end

type SuperMatrix
  diml_::Int64
  dimr_::Int64
  dim_::Int64
  opl_::Ops
  opr_::Ops
end


import Base.issymmetric
function issymmetric(x::SuperMatrix)
  return true
end


import Base.eltype
function eltype(x::SuperMatrix)
  return eltype(rand(2));
end

import Base.A_mul_B!
function A_mul_B!(v2_, x::SuperMatrix, v1_)
  v1 = reshape(v1_, x.diml_, x.dimr_);
	v2 = x.opl_.heff * v1;
	v2 = v2 + v1 * x.opr_.heff;

	v_1 = 0.5 * x.opl_.sp * v1 * x.opr_.sp;
  v_2 = 0.5 * x.opl_.sp' * v1 * x.opr_.sp';
	v2 = v2 + v_1 + v_2 + x.opl_.sz * v1 * x.opr_.sz;
  v2 = reshape(v2, x.dim_);
  v2_[:] = v2[:];
end


import Base.size
function size(x::SuperMatrix)
  return x.dim_,x.dim_
end


function ED_1D_hessenberg(latl::Int64, jz::Float64)
sp_ = [0. 1.;0. 0.];
sz_ = [0.5*sqrt(jz) 0.;0. -0.5*sqrt(jz)];
heff_ = zeros(2,2);
eye_ = eye(2,2);
op_l = Ops(eye_, sp_, sz_, heff_);
half_l = convert(Int64, latl/2);
dim = 2;
for i_l in 1:1:half_l-1
	dim = 2^(i_l + 1);
	h_i = 0.5*kron(op_l.sp, sp_');
	op_l.heff = kron(op_l.heff,eye(2,2)) + h_i + h_i' + kron(op_l.sz,sz_);
	op_l.sp = kron(op_l.eye, sp_);
	op_l.sz = kron(op_l.eye, sz_);
	op_l.eye = eye(dim, dim);
end
op_r = op_l;

# ED
ED_object = SuperMatrix(dim, dim, dim*dim, op_l, op_r);
E_gs = eigs(ED_object, nev=1, which=:SR);
return E_gs[1];
end





end # module
