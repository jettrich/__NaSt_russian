<TeXmacs|1.0.7.18>

<style|generic>

<\body>
  Solve Navier-Stokes equation

  <\eqnarray*>
    <tformat|<table|<row|<cell|<choice|<tformat|<table|<row|<cell|>|<cell|<frac|d
    ee|d t>=-u*<frac|d ee|d x>-v <frac|d ee|d y>\<noplus\>+Re*\<Delta\>
    <around*|(|ee|)>>>|<row|<cell|>|<cell|ee=\<Delta\><around*|(|ps|)>>>|<row|<cell|>|<cell|<frac|d
    ps|d y>=u<rsub|>\<nocomma\>,<frac|d ps|d x>=-v>>>>>>|<cell|>|<cell|>>>>
  </eqnarray*>

  with border conditions:

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<downarrow\>O
    y>|<cell|B3>|<cell|\<rightarrow\>O x>>|<row|<cell| input
    B4>|<cell|>|<cell|B6 output>>|<row|<cell|>|<cell|B1>|<cell|>>>>
  </eqnarray*>

  \ B4 is steam input. It weak near edges (coasts) and strong near center
  \ so it use parabola equation for speed <math|V=(u(y),0) >where<math|
  u(y)=C*y*(H-y)>. For ps used :

  <\equation*>
    \ ps<around*|(|y|)>=<big|int><rsub|0><rsup|y>u(s) d s=<frac|C H
    y<rsup|2>|2>-<frac|C H y<rsup|3>|3>+const
  </equation*>

  \ 

  B6 is same like B4

  \ B1 and B3 we use <math|u=<frac|d ps|d y>,v=-<frac|d ps|d x>
  \<Rightarrow\>ps=v*x+c=c1*x+c=u0*x+c> (u0 is const for init)

  for <math|ee> in all borders Bi \ used deffinition of function of vorticity
  \ 

  <\equation*>
    ee=<frac|d u|d y>-<frac|d v|d x>
  </equation*>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>