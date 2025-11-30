





po <- po_simul(m1 = 10, m0 = 20)

# true ate
10 - 20

mean(po$delta)

sd(po$delta)


mean(
  po$Y[po$A == 1]

) - mean(
  po$Y[po$A == 0]
)



