FDTD EM Wave Simulations in 2D
==============

A 3rd year project at KCL..

## Requirements
numpy, matplotlib, json

## Usage
See ```tutorial.py``` for full example.

#### Changing boundary properties
Change boundaries of simulation to reflective walls (default walls absorb).  
Note that calling the function without arguments will set all walls to reflective.
```python
solver.set_reflect_boundaries(up=True, down=True, left=False, right=False)
```

#### Editing material
The upper left and lower right parameters specify an iterable which contains the coordinates of the
upper left corner and the lower right corner of the rectangle.

After adding all pulses, user can create material matrix (the matrix containing permittivity and permeability values):
```python
material = solver.create_material()
```
This object has methods that allows adding different shaped materials. For example:
```python
material.set_material_rect(upper_left=(x1, y1), lower_right=(x2, y2))
material.set_material_convex(centre=(1.5, 1.5), radius=2, thickness=1, epsilon_rel=5)
```
A plot of the materials can be shown using the following method:
```python
material.plot_time()
```

## TODO

* Fix waveguide 
