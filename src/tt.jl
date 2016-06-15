module tt
global e1=1

function t()
  global e1 = 2
end

function s()
  global e1 = e1+1
end
end
