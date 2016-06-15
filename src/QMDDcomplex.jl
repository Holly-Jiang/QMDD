module QMDDcomplex

export Qmake, QMDDcos, QMDDsin, Clookup, QMDDinitCtable, QMDDcomplexInit, QMDDmakeRootsOfUnity

global COMPLEXTSIZE = 20000
global MAXRADIX = 2
global Radix = 2

global Ctable = Vector{Complex}(COMPLEXTSIZE)
global Cmag = Matrix{Float64}(COMPLEXTSIZE,COMPLEXTSIZE)
global Cangle = Matrix{Float64}(COMPLEXTSIZE,COMPLEXTSIZE)
global CTa = Vector{Int64}(MAXRADIX+1)
global Ctol = 1e-10
global Ctentries = 2

  function Qmake(a::Int64, b::Int64, c::Int64)
    return ((1.0a+1.0b)*√2.0)/c
  end

  function QMDDcos(fac::Float64, div::Float64)
    return cos(fac)*π/div
  end

  function QMDDsin(fac::Float64, div::Float64)
    return sin(fac)*π/div
  end

  function QMDDinitCtable()
    fill!(Ctable,0.0)
    Ctable[2]=1.0
    return 0
  end

  function QMDDinitComplex()
    Cmag[1]=0
    Cmag[2]=1
    QMDDinitCtable()
  end

  function angle(a::Int64)
    ca=Ctable[a]
    if imag(ca)>=0.0
      return acos(real(ca)/Cmag[a])
    else
      return 2π-acos(real(ca)/Cmag[a])
    end
  end

  function Clookup(c::Complex)
    i = findfirst(Ctable,c)
    if i != 0
      return i
    end
    
    global Ctentries = Ctentries + 1

    if Ctentries > COMPLEXTSIZE
      println("Complex table overflow")
      return 0
    end

    Ctable[Ctentries]=c
    Cmag[Ctentries]=abs(c)
    Cangle[Ctentries]=angle(Ctentries)
    return Ctentries
  end

  function Cnegative(a::Int)
    c = Ctable[a]
    c = -c
    return(Clookup(c))
  end

  function Cadd(ai::Int, bi::Int)
    if ai == 1
      return bi
    end
    if bi == 1
      return ai
    end

    a = Ctable[ai]
    b = Ctable[bi]

    r = a+b

    t=Clookup(r)
    return t
  end

  function Csub(ai::Int, bi::Int)
    if bi == 1
      return ai
    end
    
    a = Ctable[ai]
    b = Ctable[bi]

    r = a-b

    t = Clookup(r)
    return t
  end

  function Cmul(ai::Int, bi::Int)
    if ai == 2
      return bi
    end
    if bi == 2
      return ai
    end

    if ai==1 || bi==1
      return 0
    end

    a=Ctable[ai]
    b=Ctable[bi]

    r = a*b

    t=Clookup(r)
    return t
  end

  function CintMul(a::Int, bi::Int)
    r = Ctable[bi]
    r = a*r
    return(Clookup(r))
  end

  function Cdiv(ai::Int, bi::Int)
    if ai==bi
      return 1
    end

    if ai==0
      return 0
    end

    if bi==2
      return ai
    end

    a=Ctable[ai]
    b=Ctable[bi]

    r = a/b

    t=Clookup(r)
    return t
  end 

  function QMDDmakeRootsOfUnity()
    el = cos(2π/Radix) + (sin(2π/Radix))im
    CTa[1]=1
    CTa[2]=Clookup(el)
    [CTa[x]=Cmul(CTa[x-1],CTa[2]) for x=3:Radix+1]
  end

  function CAbs(a::Int)
    if a < 3
      return a
    end
    
    s = Ctable[a]

    r = Cmag[a] + 0.0im

    b = Clookup(r)

    return b
  end

  function CUnit(a::Int)
    if a < 3
      return a
    end 
    
    if Cmag[a] < 1.0
      return 0
    else
      return 1
    end
  end
end
