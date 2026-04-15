#!/usr/bin/env python

# Standard library imports
import argparse
import logging
from itertools import cycle
from pathlib import Path
from typing import List, Dict, Optional

# Third-party library imports
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import plotly
import plotly.express as px
import plotly.graph_objects as go


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Process clustering and distance data for analysis."
    )

    # Input arguments
    input_group = parser.add_argument_group("Input Arguments")
    input_group.add_argument(
        "--cluster_stat",
        type=Path,
        required=True,
        help="Path to the output TSV file containing intra- and inter-cluster distance statistics.",
    )

    input_group.add_argument(
        "--distance_to_count",
        type=Path,
        required=True,
        help="Path to the output TSV file with counts of intra- and inter-cluster distances, "
        "used for generating density plots.",
    )

    # Output arguments
    output_group = parser.add_argument_group("Output Arguments")
    output_group.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        default="plots",
    )
    return parser.parse_args(argv)


def make_box_plot_figure(df_distance_metric: pd.DataFrame) -> go.Figure:
    """
    Create a grouped box plot visualizing intra-cluster and inter-cluster distances
    for different species.

    :param df_distance_metric: DataFrame containing distance metrics for each species.
                               Expected columns:
                               - "species": Species names.
                               - "intra_cluster_q1", "intra_cluster_median", "intra_cluster_mean",
                                 "intra_cluster_q3", "intra_cluster_min", "intra_cluster_max": Metrics for intra-cluster distances.
                               - "inter_cluster_q1", "inter_cluster_median", "inter_cluster_mean",
                                 "inter_cluster_q3", "inter_cluster_min", "inter_cluster_max": Metrics for inter-cluster distances.
    :return: A Plotly Figure object with the box plot visualization.
    """
    # Define colors for intra-cluster and inter-cluster distances
    inter_dist_color = "#DC3220"  # Red for inter-cluster
    intra_dist_color = "#005AB5"  # Blue for intra-cluster

    fig = go.Figure()

    # Add intra-cluster distance box plot
    fig.add_trace(
        go.Box(
            name="Intra-cluster distance",
            x=df_distance_metric["species"],
            q1=df_distance_metric["intra_cluster_q1"],
            median=df_distance_metric["intra_cluster_median"],
            mean=df_distance_metric["intra_cluster_mean"],
            q3=df_distance_metric["intra_cluster_q3"],
            lowerfence=df_distance_metric["intra_cluster_min"],
            upperfence=df_distance_metric["intra_cluster_max"],
            marker_color=intra_dist_color,
        )
    )

    # Add inter-cluster distance box plot
    fig.add_trace(
        go.Box(
            name="Inter-cluster distance",
            x=df_distance_metric["species"],
            q1=df_distance_metric["inter_cluster_q1"],
            median=df_distance_metric["inter_cluster_median"],
            mean=df_distance_metric["inter_cluster_mean"],
            q3=df_distance_metric["inter_cluster_q3"],
            lowerfence=df_distance_metric["inter_cluster_min"],
            upperfence=df_distance_metric["inter_cluster_max"],
            marker_color=inter_dist_color,
        )
    )

    # Customize layout and appearance
    fig.update_layout(
        xaxis_title="Species",
        yaxis_title="Distances",
        template="plotly_white",
        boxmode="group",  # Group boxes for better comparison
    )

    logging.info("Box plot figure created successfully.")
    return fig


def build_density_table(
    df_count: pd.DataFrame, smoothing_parameter: float = 0.2
) -> pd.DataFrame:
    """
    Builds a density table from a count table using Gaussian kernel density estimation (KDE).

    :param df_count: DataFrame containing distances and their respective counts.
                     Expected columns:
                     - "distance": Numeric distances to be smoothed.
                     - "count": Counts associated with each distance.
    :param smoothing_parameter: Smoothing parameter for the Gaussian KDE. Smaller values lead to sharper curves.
    :return: A DataFrame containing smoothed density values and their corresponding distances.
             Columns:
             - "Distance": Normalized distances after KDE.
             - "Density": Smoothed density values.
    """
    multiplier = 10000  # Scale factor for more precise integer calculations
    df_count["distance"] = (df_count["distance"] * multiplier).astype(int)

    # Determine the maximum distance for KDE x-axis values
    max_dist = df_count["distance"].max()

    # Generate x values for KDE in the range [0, max_dist]
    x_vals = np.linspace(0, max_dist, max_dist)

    # Compute Gaussian KDE with specified smoothing parameter
    density = gaussian_kde(df_count["distance"], weights=df_count["count"])
    density.covariance_factor = lambda: smoothing_parameter
    density._compute_covariance()

    # Build density data
    density_info = {
        "Distance": x_vals / multiplier,  # Normalize distances back to original scale
        "Density": density(x_vals),
    }

    return pd.DataFrame(data=density_info)


def build_density_table_per_sp_and_type(df_count: pd.DataFrame) -> pd.DataFrame:
    """
    Builds density tables for each species and distance type using Gaussian KDE.

    :param df_count: DataFrame containing distances, species, types, and their respective counts.
                     Expected columns:
                     - "distance": Numeric distances.
                     - "count": Counts associated with each distance.
                     - "species": Species identifier.
                     - "type": Type of distance (e.g., "intra_cluster", "inter_cluster").
    :return: A combined DataFrame containing smoothed density values for all species and types.
             Columns:
             - "Distance": Normalized distances.
             - "Density": Smoothed density values.
             - "type": Distance type (e.g., "intra_cluster", "inter_cluster", "all").
             - "sp": Species identifier.
    """
    # Aggregate counts across all distance types for each species
    df_count_all = (
        df_count.groupby(["distance", "species"]).agg({"count": "sum"}).reset_index()
    )
    df_count_all["type"] = "all"

    # Combine aggregated counts with the original count table
    df_count = pd.concat([df_count, df_count_all])
    df_density_list = []

    # Build density tables for each species and distance type
    for species in df_count["species"].unique():
        for dist_type in ["intra_cluster", "inter_cluster", "all"]:
            df_count_type = df_count.loc[
                (df_count["type"] == dist_type) & (df_count["species"] == species)
            ].copy()

            # Skip if there are insufficient data points for density estimation
            if len(df_count_type["count"]) <= 50:
                continue

            df_density_type = build_density_table(df_count_type)
            df_density_type["type"] = dist_type
            df_density_type["sp"] = species

            df_density_list.append(df_density_type)

    # Combine density tables for all species and types
    if df_density_list:
        return pd.concat(df_density_list)
    else:
        logging.warning("No density tables were created due to insufficient data.")
        return pd.DataFrame(columns=["Distance", "Density", "type", "sp"])


def make_cluster_distance_density_plots_per_sp(
    df_density: pd.DataFrame,
) -> Dict[str, plotly.graph_objects.Figure]:
    """
    Creates density plots for each species, showing distributions for intra-cluster,
    inter-cluster, and overall distances.

    :param df_density: DataFrame containing density data for each species and distance category.
                       Expected columns:
                       - "sp": Species identifier.
                       - "type": Distance type ("intra_cluster", "inter_cluster", "all").
                       - "Distance": X-axis values for the density plot.
                       - "Density": Y-axis density values.
    :return: A dictionary mapping species names to their corresponding Plotly area plots.
    """
    color_discrete_map = {
        "intra_cluster": "#DC3220",
        "inter_cluster": "#005AB5",
        "all": "#6B2D5C",
    }

    sp_to_figures = {}
    df_density["Distance category"] = df_density[
        "type"
    ]  # Rename for better plot labels

    for sp in df_density["sp"].unique():
        # Filter data for the current species
        df_density_sp = df_density.loc[df_density["sp"] == sp]

        # Create area plot
        fig = px.area(
            df_density_sp,
            x="Distance",
            y="Density",
            color="Distance category",
            facet_row="Distance category",
            category_orders={
                "Distance category": ["all", "intra_cluster", "inter_cluster"]
            },
            title=f"{sp}: Distance distribution",
            color_discrete_map=color_discrete_map,
        )
        fig = fig.update_traces(line=dict(width=2), opacity=0.7)

        sp_to_figures[sp] = fig

    return sp_to_figures


def make_distance_density_plots(
    df_density: pd.DataFrame, sorted_species: List[str]
) -> Dict[str, plotly.graph_objects.Figure]:
    """
    Creates density plots for all species, showing general distance distributions.
    Generates two types of plots:
    1. Area plots with species as facets.
    2. Line plots stacked for all species.

    :param df_density: DataFrame containing density data for all species.
                       Expected columns:
                       - "sp": Species identifier.
                       - "type": Distance type ("all").
                       - "Distance": X-axis values for the density plot.
                       - "Density": Y-axis density values.
    :param sorted_species: List of species names sorted in the desired order for the plots.
    :return: A dictionary mapping output filenames to their corresponding Plotly figures.
    """
    # Define a color palette for species
    colors = (
        px.colors.qualitative.Bold
        + px.colors.qualitative.G10
        + px.colors.qualitative.Pastel
        + px.colors.qualitative.Antique
        + px.colors.qualitative.Safe
        + px.colors.qualitative.Vivid
        + px.colors.qualitative.Dark24
    )
    color_discrete_map = {sp: color for sp, color in zip(sorted_species, cycle(colors))}

    # Filter data to include only general distances
    df_density_all = df_density.loc[df_density["type"] == "all"]

    # Area plot with species as facets
    fig_sp_row = px.area(
        df_density_all,
        x="Distance",
        y="Density",
        color="sp",
        title="Distance distribution",
        color_discrete_map=color_discrete_map,
        facet_row="sp",
        category_orders={"sp": sorted_species},
    ).update_traces(line=dict(width=2), opacity=0.8)

    # Line plot for all species
    fig_line_stacked = px.line(
        df_density_all,
        x="Distance",
        y="Density",
        color="sp",
        title="Distance distribution",
        color_discrete_map=color_discrete_map,
        category_orders={"sp": sorted_species},
    )

    # Map output filenames to figures
    name_to_fig = {
        "distance_density_plot_sp_row.html": fig_sp_row,
        "distance_density_plot_lines.html": fig_line_stacked,
    }

    return name_to_fig


def main(argv: Optional[list[str]] = None) -> None:
    """
    Coordinate argument parsing and program execution.

    This function parses input arguments, processes a distance matrix and cluster composition file,
    computes cluster metrics, generates visualizations, and saves outputs.

    :param argv: List of arguments to override command-line input for testing purposes. Defaults to None.
    """
    # Configure logging
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )

    logging.info("Starting the program.")

    # Parse arguments
    args = parse_args(argv)
    output_dir = args.output_dir

    # Ensure the output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Output directory created: {output_dir}")

    # Load and process cluster statistics
    logging.info(f"Loading cluster statistics from {args.cluster_stat}")
    df_distance_metric = pd.read_csv(args.cluster_stat, sep="\t")
    df_distance_metric = df_distance_metric.sort_values(
        by=["genomes", "species"], ascending=[False, True]
    )
    sorted_species = df_distance_metric["species"].to_list()
    logging.info("Cluster statistics loaded. Species sorted by genome count.")

    # Generate and save box plot
    logging.info("Generating box plot for intra- and inter-cluster distances.")
    box_plot = make_box_plot_figure(df_distance_metric)
    box_plot_html = output_dir / "distance_box_plot.html"
    box_plot.write_html(box_plot_html)
    logging.info(f"Box plot saved to {box_plot_html}")

    # Load pairwise distances and compute density tables
    logging.info(f"Loading pairwise distance data from {args.distance_to_count}")
    df_count = pd.read_csv(args.distance_to_count, sep="\t")
    logging.info("Generating density tables for each species and distance type.")
    df_density = build_density_table_per_sp_and_type(df_count)

    # Generate and save density plots per species
    logging.info("Creating density plots for each species.")
    sp_to_figures = make_cluster_distance_density_plots_per_sp(df_density)
    for sp, figure in sp_to_figures.items():
        html_fig_path = output_dir / f"{sp}_cluster_distances_density_plot.html"
        figure.write_html(html_fig_path)
        logging.info(f"Density plot for {sp} saved to {html_fig_path}")

    # Generate and save combined density plots
    logging.info("Creating combined density plots.")
    name_to_fig = make_distance_density_plots(df_density, sorted_species)
    for name, figure in name_to_fig.items():
        html_fig_path = output_dir / name
        figure.write_html(html_fig_path)
        logging.info(f"Combined density plot saved to {html_fig_path}")

    logging.info("Program completed successfully.")


if __name__ == "__main__":
    main()
