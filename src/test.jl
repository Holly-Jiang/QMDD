using QMDDpackage

QMDDinit(1)
QMDDinitGateMatrices()

println("Hadamard: $(Hm)")

println("One QMDD node $(QMDDone)")

QMDDprint(QMDDone,10)
