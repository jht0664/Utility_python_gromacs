plot 'COM1.dnums.align.avg' u ($0/NBINS):1:2 w ye
f(x) = wp+0.5*(wr-wp)*(tanh((x-b+c)/lamda)-tanh((x-b-c)/lamda))
wr=0.85;wp=0.15;b=0.5;c=0.23;lamda=0.1
plot 'COM1.dnums.align.avg' u ($0/NBINS):1:2 w ye, f(x)
fit f(x) 'COM1.dnums.align.avg' u ($0/NBINS):1 via wr, wp, b, c, lamda
fit f(x) 'COM1.dnums.align.avg' u ($0/NBINS):1:2 via wr, wp, b, c, lamda
plot 'COM1.dnums.align.avg' u ($0/NBINS):1:2 w ye lt 4 lc 1 title "data", f(x) lw 2 lc 1 title "fit"

