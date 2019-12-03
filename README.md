2D_Solver_Test
==============

## EM Solver in 2D, official version!
Requires numpy, matplotlib

## Usage
See ```tutorial.py``` for example.

// TODO

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
* 20/11/2019 Added oscillating pulse, changed user input settings

## TODO

* Plane wave excitation
* Save images as gif functionality

## Questions

### Problems with boundary
* Definition of plane wave (in between nodes?)
* Do we apply the condition always or only while plane wave
is running
* If there's a reflected wave going from right to left, and we override it with
our plane wave going left to right, won't that reflected wave disappear?

### Measuring entropy?
1. Use reflecting walls such that energy is conserved
2. Collect all values in the matrix (calculate energy)
3. Use some method to measure the standard deviation between all points (or variance)
4. Create plot over time of this value, can possibly be entropy
