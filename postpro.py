#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# **** Plot the CSV data from result with experimental data *** #

import pandas as pd
import matplotlib.pyplot as plt
# import numpy as np
import os

# Variable to be plotted
# Change if needed. Based on the header name of the CSV files and name of legend!
plot_var_model = "phi" # phi, dpdz_inv, etc.
plot_var_exp = "phi" # phi, dpdz, etc.

# plot_y_name = "Pressure loss [Pa/m]"
plot_y_name = "Solid volume fraction"

# Optional input
exclude_mflux = [0, 5.1, 141] 

def plot_comparison(
    # model_path = "results/model/pneumatic_transport_results_1.csv",
    model_path = "results/model/pneumatic_transport_results_withPartFrict0.01.csv",
    exp_path = "results/rautiainen_1999_experiment/figure4_phi_over_Vf.csv",
    output_dir = "results/plots",
    # output_plot_filename = "model_comparison_phi.png",
    output_plot_filename = "model_comparison_phi_withPartFrict0.01.png",
    ):

    # Load model and experiment data
    model_data = pd.read_csv(model_path)
    exp_data = pd.read_csv(exp_path)

    for col in ["Vf",plot_var_model,"mflux_part"]:
        model_data[col] = pd.to_numeric(model_data[col], errors="coerce")
    model_data.dropna(subset = ["Vf",plot_var_model,"mflux_part"], inplace=True)

    model_data = model_data[~model_data["mflux_part"].isin(exclude_mflux)]

    for col in ["Vf",plot_var_exp,"mflux_part"]:
        exp_data[col] = pd.to_numeric(exp_data[col], errors="coerce")
    exp_data.dropna(subset = ["Vf",plot_var_exp,"mflux_part"], inplace=True)

    exp_data = exp_data[~exp_data["mflux_part"].isin(exclude_mflux)]


    os.makedirs(output_dir, exist_ok=True)
    outpath = os.path.join(output_dir, output_plot_filename)

    mflux_part_case = sorted(set(model_data["mflux_part"].unique()) | set(exp_data["mflux_part"].unique()))
    cmap = plt.cm.get_cmap("tab20", len(mflux_part_case))
    color_for = {mf: cmap(i) for i, mf in enumerate(mflux_part_case)}

    plt.figure(figsize=(10,8))

    for mflux, grp in model_data.sort_values(["mflux_part","Vf"]).groupby("mflux_part"):
        plt.plot(grp["Vf"], grp[plot_var_model], linestyle="-",
                 linewidth=2, color=color_for[mflux],
                 label=f"Model, {mflux:g} kg/m2/s"
                 )
        
    for mflux, grp in exp_data.sort_values(["mflux_part", "Vf"]).groupby("mflux_part"):
        plt.plot(grp["Vf"], grp[plot_var_exp], linestyle="-.", marker="x",
                 markersize=8, linewidth=2,
                 color=color_for[mflux], label=f"Exp, {mflux:g} kg/m2/s"
                )
        
    # If needed, uncomment. To annotate value of fitting coefficient, etc.
    plt.text(x=8.98, y=0.00835, s="A = 0.01", fontsize=12, color="black")
        
    plt.xlabel("Fluid velocity [m/s]")
    plt.ylabel(plot_y_name)

    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, ncol=2)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    plot_comparison(
        # model_path="results/model/pneumatic_transport_results_1.csv",
        model_path="results/model/pneumatic_transport_results_withPartFrict0.01.csv",
        exp_path="results/rautiainen_1999_experiment/figure4_phi_over_Vf.csv",
    )









