# NavierStokes
```
A simple example of flow in a tube with barrier
See description.pdf for more
```
### requirement
python2.7,numpy,scipy,pylab
### usage
```
Main class contained in basis.py
File extend.py contains additional conditions that don't tested yet
```
```
from basis import NS
ns=NS()
ns.testDefault(timeSteps=22)

it return phase field of solution after 22 steps
should look like:
```

 ![alt tag](https://raw.githubusercontent.com/valdecar/NavierStokes/master/flow.png)

### References:
```
Роуч П. Вычислительная гидродинамика
```
