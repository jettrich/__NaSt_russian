# NavierStokes
a simple example of flow in a tube with barrier

### requirement
numpy,scipy,pylab
### usage 
```
from basis import NS
ns=NS()
ns.testDefault(timeSteps=22)

it return phase field of solution after 22 steps
shold  look like:
```

 ![alt tag](https://raw.githubusercontent.com/valdecar/NavierStokes/master/flow.png)
