FDTD EM Wave Simulations in 2D
==============

A 3rd year project at KCL..

## Requirements
numpy, matplotlib

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

solver.set_material_rect(upper_left= (i3, j3), lower_right=(i4, j4), epsilon_rel, mu_rel)
```

## TODO

* Save images as gif functionality
* Save h matrix as npy file

