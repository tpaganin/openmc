import numpy as np

import openmc

def fill_3d_list(n, val):
    """
    Generates a 3D list of dimensions nxnxn filled with copies of val.

    Parameters:
    n (int): The dimension of the 3D list.
    val (any): The value to fill the 3D list with.

    Returns:
    list: A 3D list of dimensions nxnxn filled with val.
    """
    return [[[val for _ in range(n)] for _ in range(n)] for _ in range(n)]

def create_random_ray_model():
    ###############################################################################
    # Create multigroup data

        # Instantiate the energy group data
    ebins = [1e-5, 1e-1, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=ebins)
    
    # High scattering ratio means system is all scattering
    # Low means fully absorbing
    scattering_ratio = 0.5

    source_total_xs = [0.1, 0.1]
    source_mat_data = openmc.XSdata('source', groups)
    source_mat_data.order = 0
    source_mat_data.set_total([0.1, 0.2])
    #source_mat_data.set_absorption([xs * (1.0 - scattering_ratio) for xs in source_total_xs])
    source_mat_data.set_absorption([0.05, 0.1])
    source_scatter_matrix = \
            [[[0.045, 0.005],
              [0.0, 0.1]]]
    source_scatter_matrix = np.array(source_scatter_matrix)
    source_scatter_matrix = np.rollaxis(source_scatter_matrix,0,3)
    source_mat_data.set_scatter_matrix(source_scatter_matrix)
    #source_scatter_matrix = np.rollaxis(np.array([[0.045, 0.005], [0.001, 0.099]]),0,2)
    #source_scatter_matrix = source_scatter_matrix.reshape((2, 2, 1))  # Reshape to (2, 2, 1)
    #source_mat_data.set_scatter_matrix(source_scatter_matrix)

    void_total_xs = [1.0e-4, 1.0e-4]  # Example total cross-sections for two energy groups
    void_mat_data = openmc.XSdata('void', groups)
    void_mat_data.order = 0
    void_mat_data.set_total([1.0e-4, 1.0e-4])
    #void_mat_data.set_absorption([xs * (1.0 - scattering_ratio) for xs in void_total_xs])
    void_mat_data.set_absorption([0.5e-4, 0.5e-4])
    void_scatter_matrix = \
            [[[4.5e-5, 5.0e-6],
               [0.0, 0.5e-4]]]
    void_scatter_matrix = np.array(void_scatter_matrix)
    void_scatter_matrix = np.rollaxis(void_scatter_matrix,0,3)
    void_mat_data.set_scatter_matrix(void_scatter_matrix)

    #void_scatter_matrix = np.rollaxis(np.array([[4.5e-5, 5.0e-6], [1.0e-6, 9.9e-5]]),0,2) # Scatter from group 1 to 1 and 2
    #void_scatter_matrix = void_scatter_matrix.reshape((2, 2, 1))  # Reshape to (2, 2, 1)
    #void_mat_data.set_scatter_matrix(void_scatter_matrix)

    # Define XSdata for shield material
    shield_total_xs = [0.1, 0.2]  # Example total cross-sections for two energy groups
    shield_mat_data = openmc.XSdata('shield', groups)
    shield_mat_data.order = 0
    shield_mat_data.set_total([0.1, 0.2])
    #source_mat_data.set_absorption([xs * (1.0 - scattering_ratio) for xs in source_total_xs])
    shield_mat_data.set_absorption([0.05, 0.1])
    shield_scatter_matrix = \
            [[[0.045, 0.005],
              [0.0, 0.1]]]
    shield_scatter_matrix = np.array(shield_scatter_matrix)
    shield_scatter_matrix = np.rollaxis(shield_scatter_matrix,0,3)
    shield_mat_data.set_scatter_matrix(shield_scatter_matrix)
    
    
    #shield_mat_data.set_total(shield_total_xs)
    #shield_mat_data.set_absorption([xs * (1.0 - scattering_ratio) for xs in shield_total_xs])
    #shield_scatter_matrix = np.rollaxis(np.array([[0.045, 0.005], [0.001, 0.099]]),0,2)# Scatter from group 2 to 1 and 2
    #shield_scatter_matrix = shield_scatter_matrix.reshape((2, 2, 1))
    #shield_mat_data.set_scatter_matrix(shield_scatter_matrix)


    # Create MGXS Library and add XSdata
    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas([source_mat_data, void_mat_data, shield_mat_data])
    mg_cross_sections_file.export_to_hdf5()

    ###############################################################################
    # Create materials for the problem

    # Instantiate some Macroscopic Data
    source_data = openmc.Macroscopic('source')
    void_data   = openmc.Macroscopic('void')
    shield_data = openmc.Macroscopic('shield')

    # Instantiate some Materials and register the appropriate Macroscopic objects
    source_mat = openmc.Material(name='source')
    source_mat.set_density('macro', 1.0)
    source_mat.add_macroscopic(source_data)
    
    void_mat = openmc.Material(name='void')
    void_mat.set_density('macro', 1.0)
    void_mat.add_macroscopic(void_data)
    
    shield_mat = openmc.Material(name='shield')
    shield_mat.set_density('macro', 1.0)
    shield_mat.add_macroscopic(shield_data)

    # Instantiate a Materials collection and export to XML
    materials_file = openmc.Materials([source_mat, void_mat, shield_mat])
    materials_file.cross_sections = "mgxs.h5"

    ###############################################################################
    # Define problem geometry
    
    source_cell = openmc.Cell(fill=source_mat, name='infinite source region')
    void_cell = openmc.Cell(fill=void_mat, name='infinite void region')
    shield_cell = openmc.Cell(fill=shield_mat, name='infinite shield region')
    
    sub = openmc.Universe()
    sub.add_cells([source_cell])
    
    vub = openmc.Universe()
    vub.add_cells([void_cell])
    
    aub = openmc.Universe()
    aub.add_cells([shield_cell])

    # n controls the dimension of subdivision within each outer lattice element
    # E.g., n = 10 results in 1cm cubic FSRs
    n = 10
    delta = 10.0 / n
    ll = [-5.0, -5.0, -5.0]
    pitch = [delta, delta, delta]

    # create another control for a second division in one specific region such as the source.
    n2 = 10
    delta2 = 10 / n2
    ll2 = [-5.0, -5.0, -5.0] #[-2.5, -2.5, -2.5]
    pitch2 = [delta2, delta2, delta2]

    source_lattice = openmc.RectLattice()
    source_lattice.lower_left = ll2
    source_lattice.pitch = pitch2
    source_lattice.universes = fill_3d_list(n2, sub)
    
    void_lattice = openmc.RectLattice()
    void_lattice.lower_left = ll
    void_lattice.pitch = pitch
    void_lattice.universes = fill_3d_list(n, vub)
    
    shield_lattice = openmc.RectLattice()
    shield_lattice.lower_left = ll
    shield_lattice.pitch = pitch
    shield_lattice.universes = fill_3d_list(n, aub)
    
    source_lattice_cell = openmc.Cell(fill=source_lattice, name='source lattice cell')
    su = openmc.Universe()
    su.add_cells([source_lattice_cell])
    
    void_lattice_cell = openmc.Cell(fill=void_lattice, name='void lattice cell')
    vu = openmc.Universe()
    vu.add_cells([void_lattice_cell])
    
    shield_lattice_cell = openmc.Cell(fill=shield_lattice, name='shield lattice cell')
    au = openmc.Universe()
    au.add_cells([shield_lattice_cell])

    z_base = [
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [vu, vu, vu, vu, au, au],
            [vu, au, au, au, au, au],
            [vu, au, au, au, au, au],
            [vu, au, au, au, au, au],
            [vu, au, au, au, au, au],
            [su, au, au, au, au, au]
                         ]
    
    z_col = [
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au]
                         ]
    
    z_high = [
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au]
                         ]
    
    z_cap = [
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au]
                         ]

    dogleg_pattern = [
            z_base,
            z_col,
            z_col,
            z_high,
            z_cap,
            z_cap
            ]
    
    x = 60.0
    x_dim = 6

    y = 100.0
    y_dim = 10

    z = 60.0
    z_dim = 6
    
    lattice = openmc.RectLattice()
    lattice.lower_left = [0.0, 0.0, 0.0]
    lattice.pitch = [x/x_dim, y/y_dim, z/z_dim]
    lattice.universes = dogleg_pattern
    
    lattice_cell = openmc.Cell(fill=lattice)

    lattice_uni = openmc.Universe()
    lattice_uni.add_cells([lattice_cell])

    x_low  = openmc.XPlane(x0=0.0,boundary_type='reflective') 
    x_high = openmc.XPlane(x0=x,boundary_type='vacuum') 
    
    y_low  = openmc.YPlane(y0=0.0,boundary_type='reflective') 
    y_high = openmc.YPlane(y0=y,boundary_type='vacuum') 
    
    z_low  = openmc.ZPlane(z0=0.0,boundary_type='reflective') 
    z_high = openmc.ZPlane(z0=z,boundary_type='vacuum') 
    
    full_domain = openmc.Cell(fill=lattice_uni, region=+x_low & -x_high & +y_low & -y_high & +z_low & -z_high, name='full domain')

    root = openmc.Universe(name='root universe')
    root.add_cell(full_domain)

    # Create a geometry with the two cells and export to XML
    geometry = openmc.Geometry(root)

    ###############################################################################
    # Define problem settings

    # Create an initial uniform spatial source for ray integration
    lower_left = (0.0, 0.0, 0.0)
    upper_right = (x, y, z)
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
    rr_source = openmc.IndependentSource(space=uniform_dist)
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_mat], strength=2.0) 

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.batches = 100
    settings.inactive = 40
    settings.particles = 10000
    settings.solver_type = 'random ray'
    settings.run_mode = 'fixed source'
    settings.random_ray['distance_active'] = 400.0
    settings.random_ray['distance_inactive'] = 100.0
    settings.random_ray['ray_source'] = rr_source
    
    # Create an initial uniform spatial source for ray integration
    #lower_left = (-pitch, -pitch, -1)
    #upper_right = (pitch, pitch, 1)
    #uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
    #rr_source = openmc.IndependentSource(space=uniform_dist, particle="random_ray")
    
    # Create the neutron source in the bottom right of the moderator
    strengths = [0.25, 0.75] # Good - fast group appears largest (besides most thermal)
    energy_points = [1.0e-2, 1.0e1]
    energy_dist = openmc.stats.Discrete(x=energy_points,p=strengths)
    #point_source_location = openmc.stats.Point((5.0, 5.0, 5.0))
    #source = openmc.IndependentSource(energy=energy_dist, space=point_source_location, strength=1.0) # base source material
    #source 1 - base
    lower_left_src_1 = [0.0, 0.0, 0.0]
    upper_right_src_1 = [10.0, 10.0, 10.0]
    spatial_distribution_1 = openmc.stats.Box(lower_left_src_1, upper_right_src_1, only_fissionable=False)
    # source 2 - top
    #lower_left_src_9 = [50.0, 90.0, 50.0]
    #upper_right_src_9 = [60.0, 100.0, 60.0]
    #spatial_distribution_9 = openmc.stats.Box(lower_left_src_9, upper_right_src_9, only_fissionable=False)


    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_mat], strength=2.0) # works
    source_1 = openmc.IndependentSource(space=spatial_distribution_1, energy=energy_dist, strength=1.0) # works
    #source_9 = openmc.IndependentSource(space=spatial_distribution_9, energy=energy_dist, strength=1.0) # works

    #settings.source = [source, rr_source]
    settings.source = [source_1]#, source_9]
    #settings.export_to_xml()

    ###############################################################################
    # Define tallies

    # Create a mesh that will be used for tallying
    #mesh = openmc.RegularMesh()
    #mesh.dimension = (x_dim, y_dim, z_dim)
    #mesh.lower_left = (0.0, 0.0, 0.0)
    #mesh.upper_right = (x, y, z)

    # Create a mesh filter that can be used in a tally
    #mesh_filter = openmc.MeshFilter(mesh)

    # Now use the mesh filter in a tally and indicate what scores are desired
    #tally = openmc.Tally(name="Mesh tally")
    #tally.filters = [mesh_filter]
    #tally.scores = ['flux']
    #tally.estimator = 'collision'
    #tally.estimator = 'analog'

    estimator = 'tracklength'

    # Case 3A
    mesh_3A = openmc.RegularMesh()
    mesh_3A.dimension = (1, y_dim, 1)
    mesh_3A.lower_left = (0.0, 0.0, 0.0)
    mesh_3A.upper_right = (10.0, y, 10.0)
    mesh_filter_3A = openmc.MeshFilter(mesh_3A)
    
    tally_3A = openmc.Tally(name="Case 3A")
    tally_3A.filters = [mesh_filter_3A]
    tally_3A.scores = ['flux']
    tally_3A.estimator = estimator
    
    # Case 3B
    mesh_3B = openmc.RegularMesh()
    mesh_3B.dimension = (x_dim, 1, 1)
    mesh_3B.lower_left = (0.0, 50.0, 0.0)
    mesh_3B.upper_right = (x, 60.0, 10.0)
    mesh_filter_3B = openmc.MeshFilter(mesh_3B)
    
    tally_3B = openmc.Tally(name="Case 3B")
    tally_3B.filters = [mesh_filter_3B]
    tally_3B.scores = ['flux']
    tally_3B.estimator = estimator
    
    # Case 3C
    mesh_3C = openmc.RegularMesh()
    mesh_3C.dimension = (x_dim, 1, 1)
    mesh_3C.lower_left = (0.0, 90.0, 30.0)
    mesh_3C.upper_right = (x, 100.0, 40.0)
    mesh_filter_3C = openmc.MeshFilter(mesh_3C)
    
    tally_3C = openmc.Tally(name="Case 3C")
    tally_3C.filters = [mesh_filter_3C]
    tally_3C.scores = ['flux']
    tally_3C.estimator = estimator

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([tally_3A, tally_3B, tally_3C])

    ###############################################################################
    #                   Exporting to OpenMC plots.xml file
    ###############################################################################

    plot = openmc.Plot()
    plot.origin = [x/2.0, y/2.0, z/2.0]
    plot.width = [x, y, z]
    plot.pixels = [60, 100, 60]
    plot.type = 'voxel'

    # Instantiate a Plots collection and export to XML
    plot_file = openmc.Plots([plot])
    #plot_file.export_to_xml()
    
    model = openmc.model.Model()
    model.geometry = geometry
    model.materials = materials_file
    model.settings = settings
    model.xs_data = mg_cross_sections_file
    model.tallies = tallies
    model.plots = plot_file

    return model

model = create_random_ray_model()
model.export_to_model_xml()
