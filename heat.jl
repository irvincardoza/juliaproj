using Ferrite, SparseArrays
grid = generate_grid(Quadrilateral,(20,20));


dimension=2
ip=Lagrange{dimension,RefCube,1}()
qr=QuadratureRule{dimension,RefCube}(2)
cellVal= CellScalarValues(qr,ip)

dh=DofHandler(grid)

add!(dh,:u,1)
close!(dh);

K= create_sparsity_pattern(dh)
# print(K)

ch = ConstraintHandler(dh);

dirichlet_boundary= union(
    getfaceset(grid,"left"),
    getfaceset(grid,"right"),
    getfaceset(grid,"top"),
    getfaceset(grid,"bottom"),
)

dbc=Dirichlet(:u,dirichlet_boundary,(x,t)->0)
add!(ch,dbc)
close!(ch)

# Below fucntions taken from ferrite.jl library

function assemble_element!(Ke::Matrix, fe::Vector, cellvalues::CellScalarValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    # Reset to 0
    fill!(Ke, 0)
    fill!(fe, 0)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            # Add contribution to fe
            fe[i] += δu * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end

function assemble_global(cellvalues::CellScalarValues, K::SparseMatrixCSC, dh::DofHandler)
    # Allocate the element stiffness matrix and element force vector
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    # Allocate global force vector f
    f = zeros(ndofs(dh))
    # Create an assembler
    assembler = start_assemble(K, f)
    # Loop over all cels
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        # Compute element contribution
        assemble_element!(Ke, fe, cellvalues)
        # Assemble Ke and fe into K and f
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

K, f = assemble_global(cellvalues, K, dh);

apply!(K, f, ch)
u = K \ f;

vtk_grid("heat_equation", dh) do vtk
    vtk_point_data(vtk, dh, u)
end