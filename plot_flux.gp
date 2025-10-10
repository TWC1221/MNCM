 set term wxt 1 title 'Flux vs x'
 set title 'Flux vs x for different alpha'
 set xlabel 'x'
 set ylabel 'phi(x)'
 plot 'flux.dat' using 1:2 with lines title 'numerical', \
      'flux.dat' using 1:3 with lines title 'exact'
