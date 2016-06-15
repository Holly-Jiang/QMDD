module QMDDpackage

using QMDDcomplex

#/**************************************************************************************
#
#    package definitions and types
#
#    Michael Miller
#    University of Victoria
#    mmiller@cs.uvic.ca
#
#    Date: July 2008
#
#**************************************************************************************/

QMDDversion = "QMDD Package V.R1 September 2015"
global MAXSTRLEN = 11
global MAXN = 100
global MAXRADIX = 2
global MAXNEDGE = 4
global MAXNODECOUNT = 2000000
global GCLIMIT1 = 250000
global GCLIMIT_INC = 00
global MAXND = 6
global MAXDIM = 64
global NBUCKET = 8192
global HASHMASK = 8191
global CTSLOTS = 16384
global CTMASK = 16383
global COMPLEXTSIZE = 20000
global COMPLEXTMASK = 127
global TTSLOTS = 2048
global TTMASK = 2047
global MAXREFCNT = 40000
global MAXPL = 65536
global DYNREORDERLIMIT = 500
global VERBOSE = 0
global DEFINE_VARIABLES = 1 # to be overwritten by other code (CHECK/FIX!)

abstract QMDDedgeAbs

export Nm, Vm, VPm, Hm, Zm, Rm, ZEROm, Qm, Sm
export QMDDtnode, QMDDinitGateMatrices, QMDDinit
export QMDDzero, QMDDone, QMDDprint

type QMDDnode
  next::QMDDnode
  ref::UInt64
  v::Int
  renormFactor::Int64
  ident::Bool
  diag::Bool
  block::Bool
  symm::Bool
  c01::Bool
  computeSpecialMatricesFlag::Bool
  edv::Vector{QMDDedgeAbs}

  QMDDnode()=(x=new();x.ref=0;x.v=0;x.renormFactor=0;x.next=x;x.edv=Vector{QMDDedge}(MAXNEDGE);x)
end

type QMDDedge <: QMDDedgeAbs
  p::QMDDnode
  w::UInt64
  sentinel::Int64

  QMDDedge()=(x=new();x.w=0;x.p=QMDDnode();x)
end

isLast{T}(n::T) = (n == n.next)

@enum CTkind add mult kronecker reduce transpose conjugateTranspose transform c0 c1 c2 none norm createHdmSign findCnmSign findBin reduceHdm renormalize

type CTentry
  a::QMDDedge
  b::QMDDedge
  r::QMDDedge
  which::CTkind

  CTentry(a,b,r,which) = new(a,b,r,which)
  CTentry() = CTentry(QMDDedge(),QMDDedge(),QMDDedge(),none)
end

global CTable = Vector{CTentry}(CTSLOTS)

type TTentry
  n::Int64
  m::Int64
  t::Int64
  line::Vector{Int64}
  ed::QMDDedge

  TTentry(n,m,t,line,ed) = new(n,m,t,line,ed)
  TTentry()=TTentry(0,0,0,zeros(Int64,MAXN),QMDDedge())
end

global TTable = Vector{TTentry}(TTSLOTS)
global QMDDid = Vector{QMDDedge}(MAXN)

type CircuitLine
  input::ASCIIString
  output::ASCIIString
  variable::ASCIIString
  ancillary::Char
  garbage::Char
end

type QMDDrevlibDescription
  n::Int64
  ngates::Int64
  qcost::Int64
  nancillary::Int64
  ngarbage::Int64
  e::QMDDedge
  totalDC::QMDDedge
  line::Array{CircuitLine,1}
  version::ASCIIString
  inperm::ASCIIString
  outperm::ASCIIString
  ngate::Char
  cgate::Char
  tgate::Char
  fgate::Char
  pgate::Char
  vgate::Char
  kind::Array{Char,1}
  dc::Array{Char,1}
  name::Array{Char,1}
  no::Array{Char,1}
  modified::Char
end

QMDDorder=Array{Float64,1}(MAXN)

ZEROm = Matrix{Complex}(MAXRADIX,MAXRADIX)
Qm = Matrix{Complex}(MAXRADIX,MAXRADIX)
Nm = Matrix{Complex}(MAXRADIX,MAXRADIX)
Hm = Matrix{Complex}(MAXRADIX,MAXRADIX)
Sm = Matrix{Complex}(MAXRADIX,MAXRADIX)
Vm = Matrix{Complex}(MAXRADIX,MAXRADIX)
VPm = Matrix{Complex}(MAXRADIX,MAXRADIX)

function QMDDinitGateMatrices()

  v = Qmake(1,0,2) + Qmake(1,0,2)im
  vc = Qmake(1,0,2) + Qmake(-1,0,2)im

  ZEROm[1,1] = Qmake(1,0,1) + 0.0im
  ZEROm[2,2] = 0.0 + 0.0im
  ZEROm[1,2] = 0.0 + 0.0im
  ZEROm[2,1] = Qmake(1,0,1) + 0.0im

  Qm[1,1] = Qmake(1,0,1) + 0.0im
  Qm[1,2] = 0.0 + 0.0im
  Qm[2,1] = 0.0 + 0.0im
  Qm[2,2] = 0.0 + 0.0im

  Nm[1,1] = 0.0+0.0im
  Nm[2,2] = 0.0+0.0im
  Nm[1,2] = Qmake(1,0,1) + 0.0im
  Nm[2,1] = Qmake(1,0,1) + 0.0im

  Hm[1,1] = Qmake(0,1,2) + 0.0im
  Hm[1,2] = Qmake(0,1,2) + 0.0im
  Hm[2,1] = Qmake(0,1,2) + 0.0im
  Hm[2,2] = Qmake(0,-1,2) + 0.0im

  Sm[1,1] = Qmake(1,0,1) + 0.0im
  Sm[2,2] = 0.0 + Qmake(1,0,1)im
  Sm[1,2] = 0.0 + 0.0im
  Sm[2,1] = 0.0 + 0.0im

  Vm[1,1] = v
  Vm[2,2] = v
  Vm[1,2] = vc
  Vm[2,1] = vc

  VPm[1,1] = vc
  VPm[2,2] = vc
  VPm[1,2] = v
  VPm[2,1] = v

  return 0
end

if DEFINE_VARIABLES == 1
  global Radix = 2
  global Nedge = 4
  global GCswitch = 1
  global Smode = 1
  global RMmode = 0
  global MultMode = 0
  global RenormalizationNodeCount = 0
  global blockMatrixCounter = 0
  global globalComputeSpecialMatricesFlag = 1
  global dynamicReorderingTreshold = DYNREORDERLIMIT
  global largestRefCount = 0
end

function QMDDinit(verb::Int)
  if verb == 1
    println(QMDDversion)
    pp = sizeof(QMDDnode)+4*sizeof(QMDDedge)
    println("Edge size $(sizeof(QMDDedge)) bytes")
    println("Node size $(pp) bytes")
    println("Max variables $(MAXN)")
  end
    global Nedge = Radix * Radix
    global QMDDnodecount = 0
    global QMDDpeaknodecount = 0
    global Nlabel = 0

    global QMDDtnode = QMDDnode()
    QMDDtnode.ident = true
    QMDDtnode.diag = true
    QMDDtnode.block = false
    QMDDtnode.symm = true
    QMDDtnode.c01 = false
    QMDDtnode.renormFactor = 1
    QMDDtnode.computeSpecialMatricesFlag = false
    QMDDtnode.v = -1

    global QMDDzero = QMDDmakeTerminal(0.0 + 0.0im)
    global QMDDone = QMDDmakeTerminal(1.0 + 0.0im)

    fill!(TTable,TTentry())
    fill!(CTable,CTentry())
    fill!(QMDDid,QMDDedge())

    return 0
end

function QMDDmakeTerminal(c::Complex)
  ed = QMDDedge()

  ed.p = QMDDtnode
  ed.w = c
  ed.sentinel = 0

  return ed
end

function QMDDmakeNonterminal(v::Int8, edv::Vector{QMDDedge})
  ed = QMDDedge()
  reduntant = true
  ed = edv[1]
  i = 2
  while reduntant && (i<Nedge)
    reduntant = isLast(edv[i].p) || ((edv[i].w==ed.e) && (edv[i].p=e.p))
    i=i+1
  end
  if reduntant
    return edv[1]
  end

  ed = QMDDnode()
  ed.w = 1
  ed.sentinel = 0
  ed.p.v = v
  ed.p.renormFactor = 1
  ed.p.computeSpecialMatricesFlag = globalComputeSpecialMatricesFlag
  #QMDDnormalize!(ed)
  #QMDDutLookup!(ed)
  return ed
end

function QMDDprint(ed::QMDDedge, limit::Int)
  n = 0
  i = 0
  local it = ed.p
  println("top edge weight $(ed.w)")

  while !isLast(it)
    println(" $(it.ref)")
    it=it.next
  end
end

function QMDDident(x::Int, y::Int)
  ed = QMDDedge()
  f = QMDDedge()
  edv = Vector{QMDDedge}(MAXNEDGE)

  if x == 0 || y == 0
    println("ERRORE")
    return QMDDzero
  end

  if y < 0
    return QMDDone
  end

  if x == 1 && !isLast(QMDDid[y].p)
    return QMDDid[y]
  end

  f = QMDDid[y].p
  if y >= 1 && (!isLast(f))
    for i in 1..Radix
      for j in 1..Radix
        if i == j
          edv[i*Radix+j] = f
        else
          edv[i*Radix+j] = QMDDzero
        end
      end
    end
    ed = QMDDmakeNonterminal(QMDDorder[y],edv)
    QMDDid[y]=ed
    return ed
  end
  for i in 1..Radix
    for j in 1..Radix
      if i == j
        edv[i*Radix+j] = QMDDone
      else
        edv[i*Radix+j] = QMDDzero
      end
    end
  end
  ed = QMDDmakeNonterminal(QMDDorder[x],edv)
  k = x+1
  while k <= y
    for i in 1..Radix
      for j in 1..Radix
        if i == j
          edv[i*Radix+j] = ed
        else
          edv[i*Radix+j] = QMDDzero
        end
      end
    end
    ed = QMDDmakeNonterminal(QMDDorder[x],edv)
    k = k + 1
  end
  if (x==0)
    QMDDid[y]=ed
  end
  return ed
end

end
