module qcost

export gate_qcost

##*************************************************************************
## Computes gate cost based on the cost table on
## Dmitri Maslov's page http://webhome.cs.uvic.ca/~dmaslov/definitions.html
## matches the cost function in QUIVER
##

TOFFOLI_GATE = 1
FREDKIN_GATE = 2
PERES_GATE = 3
INV_PERES_GATE = 4
ANCILLARY = 0

function gate_qcost(size::Int64, n::Int64, kind::Int64)
  avail::Int64=0
  cost::BigInt=0

  if kind==PERES_GATE||kind==INV_PERES_GATE 
    return 4
  end
  avail=n-size+ANCILLARY
  if size==1 
    cost=1
  elseif size==2
    cost=1
  elseif size==3 
    cost=5
  elseif size==4
    cost=13
  elseif size==5
    if avail>=2 
      cost=26
    else 
      cost=29
    end
  elseif size==6
    if avail>=3 
      cost=38
    elseif avail>=1 
      cost=52
    else 
      cost=6
    end
  elseif size==7
    if avail>=4 
      cost=50
    elseif avail>=1 
      cost=80
    else 
      cost = 125
    end
  elseif size==8
    if avail>=5 
      cost=62
    elseif avail>=1
      cost=100
    else 
      cost=25
    end
  elseif size==9
    if avail>=6 
      cost=74
    elseif avail>=1
      cost=128
    else 
      cost=509
    end
  elseif size==10
    if avail>=7
      cost=86
    elseif avail>=1
      cost=152
    else 
      cost=1021
    end
  else 
    if avail>=size-3
      cost=12*size-34
    elseif avail>=1 
      cost=24*size-88
    else 
      cost=(BigInt(1)<<size)-3
    end
  end
  return cost
end
end
