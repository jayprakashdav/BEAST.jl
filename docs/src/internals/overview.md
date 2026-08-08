
# [Internals](@id InternalsRef)


- In general, the framework is designed such that it allows to easily add support for more kernels, finite element spaces, and excitations.
- Key are assembly routines that take in symbolic representations of the defining bilinear form. Support for block systems and finite element spaces defined in terms of direct products or tensor products of atomic spaces.


```@meta
CurrentModule = BEAST
```


## Basis

Sets of both trial and testing functions are implemented by models following the basis concept. The term basis is somewhat misleading as it is nowhere required nor enforced that these functions are linearly independent. Models implementing the Basis concept need to comply to the following semantics.


- [`numfunctions(basis)`](@ref numfunctions): number of functions in the Basis.
- [`coordtype(basis)`](): type of (the components of) the values taken on by the functions in the Basis.
- [`scalartype(d)`](@ref): the scalar field underlying the vector space the basis functions take value in.
- [`refspace(basis)`](@ref): returns the ReferenceSpace of local shape functions on which the Basis is built.
- [`assemblydata(basis)`](@ref): `assemblydata` returns an iterable collection `elements` of geometric elements and a look table `ad` for use in assembly of interaction matrices. In particular, for an index `element_idx` into `elements` and an index `local_shape_idx` in basis of local shape functions `refspace(basis)`, `ad[element_idx, local_shape_idx]` returns the iterable collection of `(global_idx, weight)` tuples such that the local shape function at `local_shape_idx` defined on the element at `element_idx` contributes to the basis function at `global_idx` with a weight of `weight`.
- [`geometry(basis)`](@ref): returns an iterable collection of Elements. The order in which these Elements are encountered corresponds to the indices used in the assembly data structure.


## Reference Space

The *reference space* concept defines an API for working with spaces of local shape functions. The main role of objects implementing this concept is to allow specialization of the functions that depend on the precise reference space used.

The functions that depend on the type and value of arguments modeling *reference space* are:

- [`numfunctions(refspace, domain)`](@ref): returns the number of shape functions on each element.

## Kernel

A kernel is a fairly simple concept that mainly exists as part of the definition of a Discrete Operator. A kernel should obey the following semantics:

In many function definitions the kernel object is referenced by `operator` or something similar. This is a misleading name as an operator definition should always be accompanied by the domain and range space.

## Discrete Operator

Informally speaking, a Discrete Operator is a concept that allows for the computation of an interaction matrix. It is a kernel together with a test and trial basis. A Discrete Operator can be passed to `assemble` and friends to compute its matrix representation.

A discrete operator is a triple `(kernel, test_basis, trial_basis)`, where `kernel` is a Kernel, and `test_basis` and `trial_basis` are Bases. In addition, the following expressions should be implemented and behave according to the correct semantics:

- [`quaddata(operator,test_refspace,trial_refspace,test_elements,trial_elements)`](@ref): create the data required for the computation of element-element interactions during assembly of discrete operator matrices.
- [`integrate!(operator,test_refspace,trial_refspace,p,test_element,q_trial_element,qd, qs, out, test_space, tptr, trial_space, bptr)`](@ref): this is the single generic function, overloaded twice over. One method, dispatching on the quadrature strategy, builds an integration strategy object `qr` describing (by its type and data fields) how to compute the interaction for the given pair of elements, using data precomputed in `qd`; the indices `p` and `q` refer to the position of the elements in the enumeration defined by `geometry(basis)` and allow fast retrieval of the relevant pre-stored data. Rather than returning `qr`, that method immediately calls the *other* method from within the same method/branch, which computes the local interaction matrix into the target buffer `zlocal`. Building and consuming `qr` in the same branch like this, instead of returning it to a separately-compiled caller, is what avoids a dynamic dispatch on `qr`'s type (which depends on the runtime geometry of the interacting elements, so is only known at runtime). Pass `action=BEAST.ReturnQRule()` to get `qr` back unevaluated instead of the default `action=BEAST.ApplyIntegrate()`. (Before BEAST 2.10 these were two separate functions, `quadrule` and `momintegrals!`; `quadrule` remains a distinct function for a few operator families outside `IntegralOperator`, such as local operators, excitations, and farfield/nearfield postprocessing.)

In the context of fast methods such as the Fast Multipole Method other algorithms on Discrete Operators will typically be defined to compute matrix vector products. These algorithms do not explicitly compute and store the interaction matrix (this would lead to unacceptable computational and memory complexity).

```@docs; canonical=false
elements
```

```@docs; canonical=false
numfunctions
scalartype
assemblydata
geometry
refspace
```

```@docs; canonical=false
quaddata
integrate!
```
