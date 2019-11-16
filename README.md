2D_Solver_Test
==============

## EM Solver in 2D, official version!
Requires numpy, matplotlib

## Usage
See ```user_test.py``` for example.

```python
from user_solver import Solver
```

Instantiate Solver object with chosen parameters.
``` python
# note s is sampling points per wavelength
solver = Solver(sigma_w, omega_0, s, stability)
```
Call the solve method for animation to begin.
```python
solver.solve()
```
### Optional:
These can be called before using the ```solve()``` method. 

#### Changing boundary properties
Change boundaries of simulation to reflective walls (default walls absorb).  
Note that calling the function without arguments will set all walls to reflective.
```python
solver.set_reflect_boundaries(up=True, down=True, left=False, right=False)
```

#### Adding materials
A rectangular reflective or refractive object can be added by the following methods.

The upper left and lower right parameters specify an iterable which contains the coordinates of the
upper left corner and the lower right corner of the rectangle.
```python
solver.add_reflect_square(upper_left= (i1, j1), lower_right=(i2, j2))

solver.add_material_square(upper_left= (i3, j3), lower_right=(i4, j4), epsilon_rel, mu_rel)
```

## Change Log


* 30/10/2019 Intiated
* 16/11/2019 Migrated; changed to OOP style, added user usage sample

## TODO

* Plane wave excitation
* Save images as gif functionality
* Oscillating pulse (should inherit Pulse class)

## Questions

### Measuring entropy?
1. Use reflecting walls such that energy is conserved
2. Collect all values in the matrix 
3. Use some method to measure the standard deviation between all points (or variance)
4. Create plot over time of this value, can possibly be entropy
