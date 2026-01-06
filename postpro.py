#################################################################
#       - 1D Solid-Gas Flow in Pneumatic Transport -            #
#                  Written for GLiP project                     #
#                           by                                  # 
#  Gabriele SIAHAAN, Maral EBADI, Ehab bin Nauman & Wajih UDDIN #
#################################################################

# **** Plot the CSV data from result with experimental data *** #

import pandas as pd
import matplotlib.pyplot as plt
import os

# Path to files or folder
# model_path = "results/model/pneumatic_transport_results_withPartFrict0.002.csv"
model_path = "results/model/pneumatic_transport_results_C1_0.1_partFrict_0.01_2.csv"
exp_path = "results/rautiainen_1999_experiment/figure4_phi_over_Vf.csv"
output_dir = "results/plots"
# output_plot_filename = "model_comparison_phi_withPartFrict0.002.png"
output_plot_filename = "model_comparison_phi_C1_1_partFrict_0.01_new.png"

# Variable to be plotted
# Change if needed. Based on the header name of the CSV files and name of legend!
plot_xVar_model = "Vf"
plot_xVar_exp = "Vf"

plot_yVar_model = "phi" # phi, dpdz_inv, etc.
plot_yVar_exp = "phi" # phi, dpdz, etc.

plot_x_name = "Fluid velocity [m/s]"
plot_y_name = "Solid volume fraction"
# plot_y_name = "Solid volume fraction" # Vsl, Vp, etc.

# Optional input
# exclude_mflux = [5.1, 10.5, 19.6, 31.8, 49.1, 59, 80.7, 141]
exclude_mflux = [0, 141]

def plot_comparison(
    model_path = model_path,
    exp_path = exp_path,
    output_dir = output_dir,
    output_plot_filename = output_plot_filename,
    ):

    # Load model and experiment data
    model_data = pd.read_csv(model_path)
    exp_data = pd.read_csv(exp_path)

    for col in [plot_xVar_model,plot_yVar_model,"mflux_part"]:
        model_data[col] = pd.to_numeric(model_data[col], errors="coerce")
    model_data.dropna(subset = [plot_xVar_model,plot_yVar_model,"mflux_part"], inplace=True)

    model_data = model_data[~model_data["mflux_part"].isin(exclude_mflux)]

    for col in [plot_xVar_exp,plot_yVar_exp,"mflux_part"]:
        exp_data[col] = pd.to_numeric(exp_data[col], errors="coerce")
    exp_data.dropna(subset = [plot_xVar_exp,plot_yVar_exp,"mflux_part"], inplace=True)

    exp_data = exp_data[~exp_data["mflux_part"].isin(exclude_mflux)]


    os.makedirs(output_dir, exist_ok=True)
    outpath = os.path.join(output_dir, output_plot_filename)

    mflux_part_case = sorted(set(model_data["mflux_part"].unique()) | set(exp_data["mflux_part"].unique()))
    cmap = plt.cm.get_cmap("tab20", len(mflux_part_case))
    color_for = {mf: cmap(i) for i, mf in enumerate(mflux_part_case)}

    plt.figure(figsize=(10,8))

    for mflux, grp in model_data.sort_values(["mflux_part",plot_xVar_model]).groupby("mflux_part"):
        plt.plot(grp[plot_xVar_model], grp[plot_yVar_model], linestyle="-",
                 linewidth=4, color=color_for[mflux],
                 label=f"Model, {mflux:g} kg/m2/s"
                 )
        
    for mflux, grp in exp_data.sort_values(["mflux_part", plot_xVar_exp]).groupby("mflux_part"):
        plt.plot(grp[plot_xVar_exp], grp[plot_yVar_exp], linestyle="-.", marker="x",
                 markersize=14, linewidth=4,
                 color=color_for[mflux], label=f"Exp, {mflux:g} kg/m2/s"
                )
        
    # If needed, uncomment. To annotate value of fitting coefficient, etc.
    # plt.text(x=10, y=2, s="C1 = 0.1", fontsize=24, color="black")
        
    plt.xlabel(plot_x_name, fontsize=22)
    plt.ylabel(plot_y_name, fontsize=22)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, labels, ncol=2, fontsize=16, loc='upper right')
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    plot_comparison(model_path, exp_path, output_dir, output_plot_filename)









