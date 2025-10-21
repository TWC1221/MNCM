 set term qt 1 noraise title 'Multigroup Diffusion Flux Profile'
 set title 'Flux vs x for three energy groups'
 set xlabel 'x'
 set ylabel 'phi(x)'
 plot 'flux.dat' using 1:2 with lines title 'g=1',\
      'flux.dat' using 1:3 with lines title 'g=2',\
      'flux.dat' using 1:4 with lines title 'g=3'
