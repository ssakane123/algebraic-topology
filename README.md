# algebraic-topology

algebraic-topology is a library to compute objects on algebraic topology, for instance, finite simplicial complex homology groups, Betti numbers, so on and so forth. 

algebraic-topology has the numpy library. 
You should install it first. 

## Usage
Here is a shot usage: 
```python
>>> from simplicial_complex import SimplicialComplex
>>> X = SimplicialComplex([[0], [1], [2], [0, 1], [0, 2], [1, 2]])
>>> X.chain()
C[0] = Z³
C[1] = Z³
>>> X.boundary_homomorphism()
d[0]: Z³ -> 0
0
d[1]: Z³ -> Z³
[[1 0 0]
 [0 1 0]
 [0 0 0]]
>>> X.cycle()
Z[0]: Z³
Z[1]: Z
>>> X.boundary()
B[0]: 0
B[1]: Z²
>>> X.homology_group()
H[0]: Z
H[1]: Z
>>> X.betti_number()
b[0]: 1
b[1]: 1
```
The notations we use above code are following: 
- `X`: a simplicial complex
- `C[q]`: the q-th module of the chain complex for `X`
- `d[q]`: the q-th boundary homomorphism of the chain complex for `X`. A matrix after the boundary homomorphism `d[i]` is the matrix representation of `d[i]`
- `Z[q]`: the kernel of the boundary homomorphism `d[i]`
- `B[q]`: the image of the boundary homomorphism `d[i]`
- `H[q]`: the q-th homology group of the simplicial complex `X`
- `b[q]`: the q-th Betti number of the simplicial complex `X`

Please note that the argument of `SimplicialComplex` have to be a simplicial complex. 