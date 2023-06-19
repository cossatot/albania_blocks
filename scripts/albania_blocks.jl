using Revise

using Oiler

using PyPlot


albania_block_file = "../block_data/albania_blocks.geojson"


daug_vels_file = "/home/itchy/research/geodesy/global_block_comps/eur_blocks/strain_data/daugostino_vels.geojson"
gsrm_vels_file = "../geod_data/gps_eur.geojson"

block_df = Oiler.IO.gis_vec_file_to_df(albania_block_file)

faults = []

@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(b->Oiler.Boundaries.boundary_to_vels(b, ee=1.0, en=1.0), 
                      non_fault_bounds)...)

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
daug_vel_df = Oiler.IO.gis_vec_file_to_df(daug_vels_file)

@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

@time daug_vels = Oiler.IO.make_vels_from_gnss_and_blocks(daug_vel_df, block_df;
    name=:site, fix="1111",
   )

gnss_vels = vcat(daug_vels,
                 gsrm_vels,
                 )

vels = vcat(#fault_vels,
            bound_vels,
            gnss_vels, 
            #geol_slip_rate_vels, 
            )

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
@info "solving"
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                                 tris=[],
                                                sparse_lhs=true,
                                                weighted=true,
                                                elastic_floor=1e-2,
                                                regularize_tris=true,
                                                tri_priors=false,
                                                tri_distance_weight=20.,
                                                predict_vels=true,
                                                check_closures=true,
                                                pred_se=false,
                                                constraint_method="kkt_sym",
                                                factorization="lu")

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)

show()

@info "writing web viewer"
Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer")
