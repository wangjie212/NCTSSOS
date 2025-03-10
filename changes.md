= TODOs
- [x] Remove unnecessary exports
- [x] Group Exports together in the NCTSSOS.jl file
- [x] Change some TODOs to NOTEs
- [x] Remove clique in MomentMethod
- [x] Remove all `DP` usage, just import with `using DynamicPolynomials: foo`

- [ ] Verify no un-necessary files in the PR

= API Redesign

```mermaid
flowchart TD
    A[Polynomial Optimization Problem] -->|_moment_relax_| B[Moment Problem]
    B -->|_apply_sparsity_| C[Moment Problems]
    C -->|_sos_dualize_| D[SOS Problem]
```

= API Naming Changes

- `sym_canon` -> `symmetric_canonalize`
- `get_ncbasis` -> `get_basis` generalized to work on commuting and non-commuting variables.

== `get_graph`
- `get_graph` -> `get_graph`
code in `nccpop.jl` and `ncupop.jl` are exactly the same
code in `trace.jl` and `complex.jl` have not been considered

= Introduced API

-  `remove_zero_degree` removes variables with degree $0$ in a monomial. It is
required for comparison of two monomials' equality.
